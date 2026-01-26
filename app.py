"""
Microbiome Dataset Discovery Dashboard

A tool for searching, scoring, and analyzing microbiome datasets
from public repositories like NCBI SRA.

Supports: NCBI SRA, ENA (European Nucleotide Archive)
Technologies: Oxford Nanopore, Illumina, PacBio
"""

import streamlit as st
import pandas as pd
import requests
import xml.etree.ElementTree as ET
from typing import Optional
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime
import json
import re


# Configuration
NCBI_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENA_API_URL = "https://www.ebi.ac.uk/ena/portal/api/search"
MAX_RESULTS = 200

# Metadata harmonization ontology - standardizing disparate metadata labels
METADATA_ONTOLOGY = {
    "age": ["age", "age_years", "host_age", "age_at_collection", "patient_age", "subject_age"],
    "sex": ["sex", "gender", "host_sex", "biological_sex", "patient_sex"],
    "bmi": ["bmi", "body_mass_index", "host_bmi"],
    "disease": ["disease", "diagnosis", "condition", "disease_status", "health_status", "phenotype"],
    "treatment": ["treatment", "medication", "drug", "therapy", "intervention"],
    "timepoint": ["timepoint", "visit", "collection_time", "sampling_time", "time_point"],
    "geographic_location": ["geo_loc_name", "country", "location", "geographic_location", "region"],
    "collection_date": ["collection_date", "sample_date", "date_collected"],
    "sample_type": ["sample_type", "tissue", "isolation_source", "body_site", "specimen"],
    "subject_id": ["subject_id", "patient_id", "participant_id", "host_subject_id", "individual"],
    "antibiotics": ["antibiotics", "antibiotic_use", "abx", "antibiotic_exposure"],
    "diet": ["diet", "dietary", "nutrition", "feeding"],
}

# Clinical metadata fields for cohort studies
CLINICAL_FIELDS = list(METADATA_ONTOLOGY.keys())

# Access classification keywords
RESTRICTED_KEYWORDS = [
    "dbgap", "controlled", "restricted", "authorized", "ega",
    "protected", "consent", "irb", "hipaa", "pii"
]


def search_sra(query: str, max_results: int = 50) -> list[str]:
    """Search NCBI SRA for datasets matching the query."""
    search_url = f"{NCBI_BASE_URL}/esearch.fcgi"
    params = {
        "db": "sra",
        "term": query,
        "retmax": max_results,
        "retmode": "json",
        "usehistory": "y"
    }

    try:
        response = requests.get(search_url, params=params, timeout=30)
        response.raise_for_status()
        data = response.json()
        return data.get("esearchresult", {}).get("idlist", [])
    except requests.RequestException as e:
        st.error(f"SRA search failed: {e}")
        return []


def fetch_sra_metadata(id_list: list[str]) -> list[dict]:
    """Fetch detailed metadata for a list of SRA IDs."""
    if not id_list:
        return []

    fetch_url = f"{NCBI_BASE_URL}/efetch.fcgi"
    params = {
        "db": "sra",
        "id": ",".join(id_list),
        "rettype": "full",
        "retmode": "xml"
    }

    try:
        response = requests.get(fetch_url, params=params, timeout=60)
        response.raise_for_status()
        return parse_sra_xml(response.text)
    except requests.RequestException as e:
        st.error(f"Metadata fetch failed: {e}")
        return []


def parse_sra_xml(xml_content: str) -> list[dict]:
    """Parse SRA XML response into structured records."""
    records = []
    try:
        root = ET.fromstring(xml_content)
    except ET.ParseError:
        return records

    for exp_package in root.findall(".//EXPERIMENT_PACKAGE"):
        record = extract_experiment_data(exp_package)
        if record:
            records.append(record)

    return records


def extract_experiment_data(exp_package) -> Optional[dict]:
    """Extract relevant fields from an EXPERIMENT_PACKAGE element."""
    record = {"source_db": "NCBI SRA"}

    # Experiment info
    experiment = exp_package.find(".//EXPERIMENT")
    if experiment is not None:
        record["accession"] = experiment.get("accession", "N/A")
        record["title"] = get_text(experiment, ".//TITLE")

        platform = experiment.find(".//PLATFORM")
        if platform is not None:
            for child in platform:
                record["platform"] = child.tag
                record["instrument"] = get_text(child, "INSTRUMENT_MODEL")
                break

    # Classify sequencing type
    platform_upper = record.get("platform", "").upper()
    instrument_upper = record.get("instrument", "").upper()

    if "OXFORD" in platform_upper or "NANOPORE" in platform_upper or "NANOPORE" in instrument_upper:
        record["seq_type"] = "Long-Read (Nanopore)"
        record["is_nanopore"] = True
    elif "PACBIO" in platform_upper:
        record["seq_type"] = "Long-Read (PacBio)"
        record["is_nanopore"] = False
    elif "ILLUMINA" in platform_upper:
        record["seq_type"] = "Short-Read (Illumina)"
        record["is_nanopore"] = False
    else:
        record["seq_type"] = "Other"
        record["is_nanopore"] = False

    # Study info
    study = exp_package.find(".//STUDY")
    if study is not None:
        record["study_accession"] = study.get("accession", "N/A")
        record["study_title"] = get_text(study, ".//STUDY_TITLE")
        record["study_abstract"] = get_text(study, ".//STUDY_ABSTRACT")

        # Extract PubMed IDs
        pubmed_ids = []
        for link in study.findall(".//STUDY_LINK//XREF_LINK"):
            db = get_text(link, "DB").lower()
            if "pubmed" in db:
                pmid = get_text(link, "ID")
                if pmid != "N/A" and pmid not in pubmed_ids:
                    pubmed_ids.append(pmid)

        record["pubmed_ids"] = ",".join(pubmed_ids) if pubmed_ids else ""

        # Extract BioProject ID
        bioproject_id = ""
        for ext_id in study.findall(".//EXTERNAL_ID"):
            namespace = ext_id.get("namespace", "")
            if namespace == "BioProject" and ext_id.text:
                bioproject_id = ext_id.text
                break

        if not bioproject_id:
            alias = study.get("alias", "")
            if alias.startswith("PRJ"):
                bioproject_id = alias

        record["bioproject_id"] = bioproject_id
        if bioproject_id:
            record["bioproject_url"] = f"https://www.ncbi.nlm.nih.gov/bioproject/{bioproject_id}"

        # Check for restricted access indicators
        abstract_lower = record.get("study_abstract", "").lower()
        title_lower = record.get("study_title", "").lower()
        is_restricted = any(kw in abstract_lower or kw in title_lower for kw in RESTRICTED_KEYWORDS)
        record["access_type"] = "Restricted (dbGaP/Controlled)" if is_restricted else "Public"

    # Sample info
    sample = exp_package.find(".//SAMPLE")
    if sample is not None:
        record["sample_accession"] = sample.get("accession", "N/A")
        record["organism"] = get_text(sample, ".//SCIENTIFIC_NAME")

        # Extract and harmonize sample attributes
        attributes = sample.findall(".//SAMPLE_ATTRIBUTE")
        sample_attrs = {}
        for attr in attributes:
            tag = get_text(attr, "TAG").lower().strip()
            value = get_text(attr, "VALUE")
            sample_attrs[tag] = value

        record["raw_metadata_fields"] = len(attributes)
        record["sample_attributes"] = sample_attrs

        # Harmonize metadata using ontology
        harmonized = harmonize_metadata(sample_attrs)
        record["harmonized_fields"] = harmonized
        record["harmonized_count"] = len([v for v in harmonized.values() if v])

        # Check for fecal/stool sample
        sample_text = " ".join(sample_attrs.values()).lower()
        is_fecal = any(term in sample_text for term in ["fecal", "feces", "stool", "gut", "intestin"])
        record["is_fecal_sample"] = is_fecal

    # Run info
    runs = exp_package.findall(".//RUN")
    total_spots = 0
    total_bases = 0
    run_accessions = []

    for run in runs:
        run_acc = run.get("accession", "")
        if run_acc:
            run_accessions.append(run_acc)
        try:
            total_spots += int(run.get("total_spots", "0"))
            total_bases += int(run.get("total_bases", "0"))
        except ValueError:
            pass

    record["run_count"] = len(runs)
    record["total_spots"] = total_spots
    record["total_bases"] = total_bases
    record["avg_read_length"] = total_bases / total_spots if total_spots > 0 else 0
    record["total_gb"] = total_bases / 1e9
    record["run_accession"] = run_accessions[0] if run_accessions else ""

    if record["run_accession"]:
        record["run_url"] = f"https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc={record['run_accession']}&display=metadata"

    if record.get("accession") and record["accession"] != "N/A":
        record["sra_url"] = f"https://www.ncbi.nlm.nih.gov/sra/{record['accession']}"

    # Library info
    library = exp_package.find(".//LIBRARY_DESCRIPTOR")
    if library is not None:
        record["library_strategy"] = get_text(library, "LIBRARY_STRATEGY")
        record["library_source"] = get_text(library, "LIBRARY_SOURCE")

    return record


def get_text(element, path: str) -> str:
    """Safely extract text from an XML element."""
    if element is None:
        return "N/A"
    found = element.find(path)
    return found.text.strip() if found is not None and found.text else "N/A"


def harmonize_metadata(raw_attrs: dict) -> dict:
    """
    Harmonize raw metadata into standardized ontology fields.
    This solves the "metadata chaos" problem by mapping disparate labels.
    """
    harmonized = {}

    for standard_field, aliases in METADATA_ONTOLOGY.items():
        value = None
        for alias in aliases:
            # Check for exact match or partial match in attribute names
            for raw_key, raw_value in raw_attrs.items():
                if alias in raw_key.lower() and raw_value and raw_value != "N/A":
                    value = raw_value
                    break
            if value:
                break
        harmonized[standard_field] = value

    return harmonized


def calculate_quality_score(record: dict, scoring_mode: str = "auto") -> dict:
    """
    Algorithmic ranking system that rates datasets based on quality metrics.
    This scoring decides what enters the curated database.
    """
    scores = {}

    # Determine scoring profile
    is_long_read = record.get("is_nanopore", False) or "PACBIO" in record.get("platform", "").upper()
    use_long_read_scoring = scoring_mode == "nanopore" or (scoring_mode == "auto" and is_long_read)
    record["scoring_mode"] = "Long-Read" if use_long_read_scoring else "Short-Read"

    if use_long_read_scoring:
        # Throughput score (0-25)
        total_gb = record.get("total_gb", 0)
        if total_gb >= 10:
            scores["depth_score"] = 25
        elif total_gb >= 5:
            scores["depth_score"] = 20
        elif total_gb >= 1:
            scores["depth_score"] = 15
        elif total_gb >= 0.1:
            scores["depth_score"] = 10
        else:
            scores["depth_score"] = 5

        # Read length score (0-25) - critical for Nanopore
        avg_length = record.get("avg_read_length", 0)
        if avg_length >= 10000:
            scores["length_score"] = 25
        elif avg_length >= 5000:
            scores["length_score"] = 20
        elif avg_length >= 1000:
            scores["length_score"] = 15
        else:
            scores["length_score"] = 10
    else:
        # Short-read scoring
        spots = record.get("total_spots", 0)
        if spots >= 100000:
            scores["depth_score"] = 25
        elif spots >= 10000:
            scores["depth_score"] = 20
        elif spots >= 1000:
            scores["depth_score"] = 15
        else:
            scores["depth_score"] = 10

        avg_length = record.get("avg_read_length", 0)
        if avg_length >= 250:
            scores["length_score"] = 20
        elif avg_length >= 150:
            scores["length_score"] = 15
        else:
            scores["length_score"] = 10

    # Metadata quality score (0-25) - harmonized fields matter more
    harmonized_count = record.get("harmonized_count", 0)
    if harmonized_count >= 8:
        scores["metadata_score"] = 25
    elif harmonized_count >= 5:
        scores["metadata_score"] = 20
    elif harmonized_count >= 3:
        scores["metadata_score"] = 15
    else:
        scores["metadata_score"] = 10

    # Clinical relevance score (0-15)
    harmonized = record.get("harmonized_fields", {})
    clinical_present = sum(1 for f in ["disease", "treatment", "age", "sex", "bmi"]
                          if harmonized.get(f))
    if clinical_present >= 4:
        scores["clinical_score"] = 15
    elif clinical_present >= 2:
        scores["clinical_score"] = 10
    else:
        scores["clinical_score"] = 5

    # Sample relevance (fecal bonus for gut microbiome)
    scores["sample_score"] = 10 if record.get("is_fecal_sample") else 5

    # Publication bonus
    scores["publication_score"] = 5 if record.get("pubmed_ids") else 0

    # Calculate total and grade
    total_score = sum(scores.values())
    scores["total_score"] = total_score

    if total_score >= 85:
        scores["quality_grade"] = "A"
        scores["ingestion_decision"] = "AUTO-INGEST"
    elif total_score >= 70:
        scores["quality_grade"] = "B"
        scores["ingestion_decision"] = "AUTO-INGEST"
    elif total_score >= 55:
        scores["quality_grade"] = "C"
        scores["ingestion_decision"] = "REVIEW"
    else:
        scores["quality_grade"] = "D"
        scores["ingestion_decision"] = "REJECT"

    record.update(scores)
    return record


def create_results_dataframe(records: list[dict]) -> pd.DataFrame:
    """Create a pandas DataFrame from the scored records."""
    if not records:
        return pd.DataFrame()

    display_columns = [
        "run_accession", "run_url", "accession", "sra_url",
        "bioproject_id", "bioproject_url", "pubmed_ids",
        "title", "organism", "platform", "instrument", "seq_type",
        "library_strategy", "total_spots", "total_gb", "avg_read_length",
        "access_type", "ingestion_decision",
        "raw_metadata_fields", "harmonized_count", "is_fecal_sample",
        "total_score", "quality_grade",
        "depth_score", "length_score", "metadata_score", "clinical_score",
        "source_db"
    ]

    df = pd.DataFrame(records)
    available_columns = [col for col in display_columns if col in df.columns]
    df = df[available_columns]

    if "total_score" in df.columns:
        df = df.sort_values("total_score", ascending=False)

    return df


def render_summary_header(df: pd.DataFrame):
    """Render the summary header with live statistics."""
    st.markdown("## Discovery Summary")

    if df.empty:
        col1, col2, col3, col4, col5 = st.columns(5)
        with col1:
            st.metric("Datasets Identified", 0)
        with col2:
            st.metric("Ready for Ingestion", 0)
        with col3:
            st.metric("Pending Review", 0)
        with col4:
            st.metric("Public Access", 0)
        with col5:
            st.metric("Restricted Access", 0)
        return

    col1, col2, col3, col4, col5 = st.columns(5)

    with col1:
        st.metric("Datasets Identified", len(df))

    with col2:
        if "ingestion_decision" in df.columns:
            auto_ingest = len(df[df["ingestion_decision"] == "AUTO-INGEST"])
            st.metric("Ready for Ingestion", auto_ingest,
                     delta=f"{auto_ingest/len(df)*100:.0f}%" if len(df) > 0 else "0%")

    with col3:
        if "ingestion_decision" in df.columns:
            review = len(df[df["ingestion_decision"] == "REVIEW"])
            st.metric("Pending Review", review)

    with col4:
        if "access_type" in df.columns:
            public = len(df[df["access_type"] == "Public"])
            st.metric("Public Access", public)

    with col5:
        if "access_type" in df.columns:
            restricted = len(df[df["access_type"].str.contains("Restricted", na=False)])
            st.metric("Restricted Access", restricted)


def render_quality_charts(df: pd.DataFrame):
    """Render quality distribution charts."""
    if df.empty:
        st.warning("No data available for charts.")
        return

    col1, col2 = st.columns(2)

    with col1:
        if "ingestion_decision" in df.columns:
            decision_counts = df["ingestion_decision"].value_counts()
            fig = px.pie(
                values=decision_counts.values,
                names=decision_counts.index,
                title="Ingestion Decision Distribution",
                color=decision_counts.index,
                color_discrete_map={
                    "AUTO-INGEST": "#2ecc71",
                    "REVIEW": "#f39c12",
                    "REJECT": "#e74c3c"
                }
            )
            st.plotly_chart(fig, use_container_width=True)

    with col2:
        if "access_type" in df.columns:
            access_counts = df["access_type"].value_counts()
            fig = px.pie(
                values=access_counts.values,
                names=access_counts.index,
                title="Public vs Restricted Access",
                color=access_counts.index,
                color_discrete_map={
                    "Public": "#3498db",
                    "Restricted (dbGaP/Controlled)": "#e74c3c"
                }
            )
            st.plotly_chart(fig, use_container_width=True)

    # Score distribution
    col1, col2 = st.columns(2)

    with col1:
        if "quality_grade" in df.columns:
            grade_counts = df["quality_grade"].value_counts().sort_index()
            fig = px.bar(
                x=grade_counts.index,
                y=grade_counts.values,
                title="Quality Grade Distribution",
                labels={"x": "Grade", "y": "Count"},
                color=grade_counts.index,
                color_discrete_map={
                    "A": "#2ecc71", "B": "#3498db",
                    "C": "#f39c12", "D": "#e74c3c"
                }
            )
            fig.update_layout(showlegend=False)
            st.plotly_chart(fig, use_container_width=True)

    with col2:
        score_cols = ["depth_score", "length_score", "metadata_score", "clinical_score"]
        available = [c for c in score_cols if c in df.columns]
        if available:
            avg_scores = df[available].mean()
            fig = px.bar(
                x=["Depth", "Read Length", "Metadata", "Clinical"],
                y=avg_scores.values,
                title="Average Score Components",
                labels={"x": "Component", "y": "Avg Score"}
            )
            st.plotly_chart(fig, use_container_width=True)

    # Technology and sample breakdown
    col1, col2 = st.columns(2)

    with col1:
        if "seq_type" in df.columns:
            tech_counts = df["seq_type"].value_counts()
            fig = px.pie(
                values=tech_counts.values,
                names=tech_counts.index,
                title="Sequencing Technology"
            )
            st.plotly_chart(fig, use_container_width=True)

    with col2:
        if "is_fecal_sample" in df.columns:
            fecal_counts = df["is_fecal_sample"].value_counts()
            fig = px.pie(
                values=fecal_counts.values,
                names=["Fecal/Gut" if x else "Other" for x in fecal_counts.index],
                title="Sample Type (Gut Relevance)"
            )
            st.plotly_chart(fig, use_container_width=True)


def render_metadata_harmonization(df: pd.DataFrame, records: list[dict]):
    """Show metadata harmonization results."""
    st.subheader("Metadata Harmonization Status")
    st.markdown(
        "Raw metadata labels are mapped to a standardized ontology. "
        "This solves the 'metadata chaos' problem for ML pipelines."
    )

    # Show ontology mapping
    with st.expander("View Harmonization Ontology"):
        ontology_df = []
        for standard, aliases in METADATA_ONTOLOGY.items():
            ontology_df.append({
                "Standard Field": standard,
                "Mapped Aliases": ", ".join(aliases[:5]) + ("..." if len(aliases) > 5 else "")
            })
        st.dataframe(pd.DataFrame(ontology_df), use_container_width=True)

    # Show harmonization coverage
    if records:
        coverage = {field: 0 for field in METADATA_ONTOLOGY.keys()}
        for rec in records:
            harmonized = rec.get("harmonized_fields", {})
            for field, value in harmonized.items():
                if value:
                    coverage[field] += 1

        coverage_df = pd.DataFrame([
            {"Field": k, "Datasets with Value": v, "Coverage %": f"{v/len(records)*100:.1f}%"}
            for k, v in coverage.items()
        ]).sort_values("Datasets with Value", ascending=False)

        st.dataframe(coverage_df, use_container_width=True)


def main():
    """Main application entry point."""
    st.set_page_config(
        page_title="Microbiome Dataset Discovery",
        page_icon=None,
        layout="wide"
    )

    st.title("Microbiome Dataset Discovery Dashboard")
    st.markdown(
        "Search, score, and analyze microbiome datasets from NCBI SRA. "
        "Supports quality scoring, metadata harmonization, and batch export."
    )

    # Sidebar
    with st.sidebar:
        st.header("Data Discovery Agent")

        st.markdown("**Source Database**")
        source_db = st.selectbox(
            "Repository",
            ["NCBI SRA", "NCBI SRA + ENA (coming soon)"],
            help="Select data source. SRA is the primary archive for sequence data."
        )

        st.markdown("---")
        st.markdown("**Search Query**")

        # Preset queries for common use cases
        preset = st.selectbox(
            "Preset Queries",
            [
                "Custom",
                "Nanopore Fecal Studies",
                "Human Gut Microbiome",
                "Clinical Stool Studies",
                "16S Amplicon (Any Platform)",
                "Shotgun Metagenomics"
            ]
        )

        preset_queries = {
            "Custom": "",
            "Nanopore Fecal Studies": "fecal[All Fields] AND Oxford Nanopore[Platform]",
            "Human Gut Microbiome": "gut microbiome[All Fields] AND human[Organism]",
            "Clinical Stool Studies": "stool[All Fields] AND clinical[All Fields] AND microbiome[All Fields]",
            "16S Amplicon (Any Platform)": "16S[All Fields] AND amplicon[All Fields] AND fecal[All Fields]",
            "Shotgun Metagenomics": "metagenome[All Fields] AND WGS[Strategy] AND gut[All Fields]"
        }

        default_query = preset_queries.get(preset, "")
        search_term = st.text_area(
            "Query",
            value=default_query if default_query else "fecal[All Fields] AND Oxford Nanopore[Platform]",
            height=80,
            help="NCBI Entrez query syntax"
        )

        max_results = st.slider("Max Results", 10, MAX_RESULTS, 50, 10)

        st.markdown("---")
        st.markdown("**Scoring Profile**")
        scoring_mode = st.radio(
            "Optimize for",
            ["auto", "nanopore", "illumina"],
            format_func=lambda x: {
                "auto": "Auto-detect",
                "nanopore": "Long-Read (Nanopore)",
                "illumina": "Short-Read (Illumina)"
            }[x]
        )

        search_button = st.button("Launch Discovery Agent", type="primary", use_container_width=True)

        st.markdown("---")
        st.markdown("**About**")
        st.caption(
            "This dashboard demonstrates automated data discovery, "
            "quality scoring, and metadata harmonization for microbiome datasets. "
            "Data sourced from NCBI SRA via E-utilities API."
        )

    # Initialize session state
    if "results_df" not in st.session_state:
        st.session_state.results_df = pd.DataFrame()
    if "records" not in st.session_state:
        st.session_state.records = []

    # Execute search
    if search_button:
        with st.spinner("Discovery Agent: Searching NCBI SRA..."):
            id_list = search_sra(search_term, max_results)

        if id_list:
            st.info(f"Agent identified {len(id_list)} potential datasets. Fetching metadata...")

            with st.spinner("Agent: Extracting and harmonizing metadata..."):
                records = fetch_sra_metadata(id_list)
                scored_records = [calculate_quality_score(r, scoring_mode) for r in records]
                st.session_state.records = scored_records
                st.session_state.results_df = create_results_dataframe(scored_records)

            st.success(f"Discovery complete. {len(st.session_state.results_df)} datasets processed and scored.")
        else:
            st.warning("No datasets found. Try adjusting your query.")

    # Display results
    df = st.session_state.results_df
    records = st.session_state.records

    # Show summary header
    render_summary_header(df)

    if not df.empty:
        st.markdown("---")

        # Tabs
        tab1, tab2, tab3, tab4 = st.tabs([
            "Dataset Queue", "Quality Analytics", "Metadata Harmonization", "Nanopore Focus"
        ])

        with tab1:
            st.subheader("Ingestion Queue")
            st.caption("Data Source: NCBI Sequence Read Archive (SRA) - https://www.ncbi.nlm.nih.gov/sra")

            # Filters
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                if "ingestion_decision" in df.columns:
                    decision_filter = st.multiselect(
                        "Ingestion Status",
                        options=df["ingestion_decision"].unique().tolist(),
                        default=df["ingestion_decision"].unique().tolist()
                    )
            with col2:
                if "access_type" in df.columns:
                    access_filter = st.multiselect(
                        "Access Type",
                        options=df["access_type"].unique().tolist(),
                        default=df["access_type"].unique().tolist()
                    )
            with col3:
                if "seq_type" in df.columns:
                    seq_filter = st.multiselect(
                        "Technology",
                        options=df["seq_type"].dropna().unique().tolist(),
                        default=df["seq_type"].dropna().unique().tolist()
                    )
            with col4:
                fecal_only = st.checkbox("Fecal/Gut samples only", value=False)

            # Apply filters
            filtered_df = df.copy()
            if "ingestion_decision" in df.columns and decision_filter:
                filtered_df = filtered_df[filtered_df["ingestion_decision"].isin(decision_filter)]
            if "access_type" in df.columns and access_filter:
                filtered_df = filtered_df[filtered_df["access_type"].isin(access_filter)]
            if "seq_type" in df.columns and seq_filter:
                filtered_df = filtered_df[filtered_df["seq_type"].isin(seq_filter)]
            if fecal_only and "is_fecal_sample" in df.columns:
                filtered_df = filtered_df[filtered_df["is_fecal_sample"] == True]

            # Display table
            st.dataframe(
                filtered_df,
                use_container_width=True,
                height=400,
                column_config={
                    "run_accession": st.column_config.TextColumn("Run ID", help="Unique SRR accession"),
                    "run_url": st.column_config.LinkColumn("Run Link", display_text="View Run"),
                    "sra_url": st.column_config.LinkColumn("Experiment", display_text="View Exp"),
                    "bioproject_url": st.column_config.LinkColumn("BioProject", display_text="View Project"),
                    "accession": st.column_config.TextColumn("Exp ID"),
                    "total_gb": st.column_config.NumberColumn("Gb", format="%.2f"),
                    "avg_read_length": st.column_config.NumberColumn("Read Len", format="%.0f"),
                    "ingestion_decision": st.column_config.TextColumn("Decision"),
                    "quality_grade": st.column_config.TextColumn("Grade"),
                }
            )

            # Export
            col1, col2 = st.columns(2)
            with col1:
                csv_data = filtered_df.to_csv(index=False)
                st.download_button(
                    "Download Full CSV",
                    csv_data,
                    "microbiome_discovery_results.csv",
                    "text/csv"
                )
            with col2:
                ingest_df = filtered_df[filtered_df["ingestion_decision"] == "AUTO-INGEST"] if "ingestion_decision" in filtered_df.columns else filtered_df
                if not ingest_df.empty:
                    st.download_button(
                        "Export Ingestion Queue",
                        ingest_df.to_csv(index=False),
                        "ingestion_queue.csv",
                        "text/csv"
                    )

            # Quick links
            st.markdown("### Source Links")
            for idx, row in filtered_df.head(10).iterrows():
                run = row.get("run_accession", "")
                run_url = row.get("run_url", "")
                bioproject = row.get("bioproject_id", "")
                bioproject_url = row.get("bioproject_url", "")
                pubmed = row.get("pubmed_ids", "")
                decision = row.get("ingestion_decision", "")
                title = str(row.get("title", ""))[:50]

                links = []
                if run_url:
                    links.append(f"[{run}]({run_url})")
                if bioproject_url:
                    links.append(f"[{bioproject}]({bioproject_url})")
                if pubmed:
                    for pmid in str(pubmed).split(","):
                        pmid = pmid.strip()
                        if pmid:
                            links.append(f"[PMID:{pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)")

                status = "✓" if decision == "AUTO-INGEST" else "?" if decision == "REVIEW" else "✗"
                if links:
                    st.markdown(f"{status} **{title}...** | {' | '.join(links)}")

        with tab2:
            st.subheader("Quality Analytics")
            render_quality_charts(df)

        with tab3:
            render_metadata_harmonization(df, records)

        with tab4:
            st.subheader("Long-Read / Nanopore Analysis")
            long_read_df = df[df["seq_type"].str.contains("Long-Read", na=False)] if "seq_type" in df.columns else pd.DataFrame()

            if not long_read_df.empty:
                st.markdown(f"**{len(long_read_df)} long-read datasets identified**")

                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.metric("Median Read Length", f"{long_read_df['avg_read_length'].median():,.0f} bp")
                with col2:
                    st.metric("Median Throughput", f"{long_read_df['total_gb'].median():.2f} Gb")
                with col3:
                    st.metric("Max Read Length", f"{long_read_df['avg_read_length'].max():,.0f} bp")
                with col4:
                    auto_ingest = len(long_read_df[long_read_df["ingestion_decision"] == "AUTO-INGEST"]) if "ingestion_decision" in long_read_df.columns else 0
                    st.metric("Auto-Ingest Ready", auto_ingest)

                # Scatter plot
                if "avg_read_length" in long_read_df.columns and "total_gb" in long_read_df.columns:
                    fig = px.scatter(
                        long_read_df,
                        x="avg_read_length",
                        y="total_gb",
                        color="quality_grade" if "quality_grade" in long_read_df.columns else None,
                        hover_data=["run_accession", "instrument"],
                        title="Read Length vs Throughput",
                        labels={"avg_read_length": "Avg Read Length (bp)", "total_gb": "Throughput (Gb)"},
                        color_discrete_map={"A": "#2ecc71", "B": "#3498db", "C": "#f39c12", "D": "#e74c3c"}
                    )
                    st.plotly_chart(fig, use_container_width=True)

                st.dataframe(long_read_df, use_container_width=True, height=300)
            else:
                st.info("No long-read datasets in current results. Try: Oxford Nanopore[Platform]")

    else:
        # Show welcome message
        st.markdown("---")
        st.markdown("""
        ### Getting Started

        1. **Select a preset query** or write a custom NCBI Entrez query
        2. **Launch the Discovery Agent** to search global repositories
        3. **Review scored datasets** with auto-ingest recommendations
        4. **Export the ingestion queue** for downstream processing

        ### Scoring Criteria

        | Component | Max Points | Description |
        |-----------|------------|-------------|
        | Sequencing Depth | 25 | Read count (Illumina) or throughput in Gb (Nanopore) |
        | Read Length | 25 | Longer reads = better taxonomic resolution |
        | Metadata Quality | 25 | Harmonized clinical/sample metadata |
        | Clinical Relevance | 15 | Disease, treatment, demographics |
        | Sample Type | 10 | Fecal/gut samples get bonus |
        | Publication | 5 | Linked PubMed ID |

        **Grades:** A (85+) = Auto-Ingest | B (70+) = Auto-Ingest | C (55+) = Review | D (<55) = Reject
        """)

    # Footer
    st.markdown("---")
    st.caption(
        "Data Discovery Engine | Sources: NCBI SRA (ncbi.nlm.nih.gov/sra) | "
        "All accessions and links verified via E-utilities API"
    )


if __name__ == "__main__":
    main()
