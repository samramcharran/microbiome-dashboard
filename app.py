"""
Microbiome Dataset Discovery Dashboard

Automated discovery and curation for microbiome research.
Searches public repositories, scores data quality, and prepares
datasets for downstream analysis pipelines.
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


# Configuration
NCBI_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MAX_RESULTS = 200

# Metadata harmonization ontology
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

# Access classification keywords
RESTRICTED_KEYWORDS = [
    "dbgap", "controlled", "restricted", "authorized", "ega",
    "protected", "consent", "irb", "hipaa", "pii"
]

# Disease categories (gut-brain axis focus)
DISEASE_CATEGORIES = {
    "Gut-Brain Axis": ["depression", "anxiety", "stress", "mood", "psychiatric", "mental",
                       "neurological", "parkinson", "alzheimer", "autism", "asd", "adhd",
                       "cognitive", "brain", "nervous system", "psycho"],
    "Pain Conditions": ["pain", "fibromyalgia", "migraine", "headache", "nociception",
                        "chronic pain", "neuropathic"],
    "GI Disorders": ["ibs", "ibd", "crohn", "colitis", "constipation", "diarrhea",
                     "gerd", "reflux", "dyspepsia", "functional gi"],
    "Metabolic": ["obesity", "diabetes", "metabolic syndrome", "insulin", "glucose",
                  "weight", "bmi", "overweight"],
    "Immune/Inflammatory": ["inflammation", "autoimmune", "allergy", "atopic", "asthma",
                            "eczema", "rheumatoid", "lupus"],
    "Infectious": ["infection", "pathogen", "cdiff", "c. difficile", "sepsis"],
    "Healthy/Control": ["healthy", "control", "normal", "reference"]
}

# Role-based preset queries
ROLE_PRESETS = {
    "Leadership": {
        "Strategic Overview": "(gut-brain[All Fields] OR depression[All Fields] OR anxiety[All Fields] OR pain[All Fields]) AND fecal[All Fields] AND microbiome[All Fields]",
        "High-Value Public Data": "fecal[All Fields] AND microbiome[All Fields] AND human[Organism]",
    },
    "Wet Lab": {
        "Cultivation Candidates (Long-Read)": "fecal[All Fields] AND Oxford Nanopore[Platform]",
        "Strain-Level Resolution": "(shotgun[All Fields] OR metagenome[All Fields]) AND fecal[All Fields] AND human[Organism]",
    },
    "Computational": {
        "ML Training Data": "gut microbiome[All Fields] AND human[Organism] AND clinical[All Fields]",
        "Large Cohorts": "stool[All Fields] AND cohort[All Fields] AND microbiome[All Fields]",
    },
    "BD / Partnerships": {
        "Collaboration Targets": "(clinical trial[All Fields] OR randomized[All Fields]) AND gut[All Fields] AND microbiome[All Fields]",
        "Probiotic Research": "probiotic[All Fields] AND (gut[All Fields] OR fecal[All Fields])",
    },
    "Disease Focus": {
        "Depression Studies": "(depression[All Fields] OR depressive[All Fields]) AND (gut[All Fields] OR fecal[All Fields]) AND microbiome[All Fields]",
        "Anxiety Studies": "anxiety[All Fields] AND (gut[All Fields] OR fecal[All Fields]) AND microbiome[All Fields]",
        "Pain Research": "(pain[All Fields] OR fibromyalgia[All Fields] OR nociception[All Fields]) AND microbiome[All Fields]",
        "IBS / IBD": "(IBS[All Fields] OR IBD[All Fields] OR \"irritable bowel\"[All Fields]) AND microbiome[All Fields]",
    }
}


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
        st.error(f"Search failed: {e}")
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

        # Check for restricted access
        abstract_lower = record.get("study_abstract", "").lower()
        title_lower = record.get("study_title", "").lower()
        is_restricted = any(kw in abstract_lower or kw in title_lower for kw in RESTRICTED_KEYWORDS)
        record["access_type"] = "Private (Requires Access)" if is_restricted else "Public"

    # Sample info
    sample = exp_package.find(".//SAMPLE")
    if sample is not None:
        record["sample_accession"] = sample.get("accession", "N/A")
        record["organism"] = get_text(sample, ".//SCIENTIFIC_NAME")

        attributes = sample.findall(".//SAMPLE_ATTRIBUTE")
        sample_attrs = {}
        for attr in attributes:
            tag = get_text(attr, "TAG").lower().strip()
            value = get_text(attr, "VALUE")
            sample_attrs[tag] = value

        record["raw_metadata_fields"] = len(attributes)
        record["sample_attributes"] = sample_attrs

        harmonized = harmonize_metadata(sample_attrs)
        record["harmonized_fields"] = harmonized
        record["harmonized_count"] = len([v for v in harmonized.values() if v])

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

    # Disease classification
    searchable_text = " ".join([
        record.get("title", ""),
        record.get("study_title", ""),
        record.get("study_abstract", "")
    ]).lower()

    detected_diseases = []
    for category, keywords in DISEASE_CATEGORIES.items():
        if any(kw in searchable_text for kw in keywords):
            detected_diseases.append(category)
    record["disease_categories"] = detected_diseases
    record["disease_category"] = detected_diseases[0] if detected_diseases else "Unclassified"
    record["gut_brain_relevant"] = "Gut-Brain Axis" in detected_diseases or "Pain Conditions" in detected_diseases

    return record


def get_text(element, path: str) -> str:
    """Safely extract text from an XML element."""
    if element is None:
        return "N/A"
    found = element.find(path)
    return found.text.strip() if found is not None and found.text else "N/A"


def harmonize_metadata(raw_attrs: dict) -> dict:
    """Harmonize raw metadata into standardized fields."""
    harmonized = {}

    for standard_field, aliases in METADATA_ONTOLOGY.items():
        value = None
        for alias in aliases:
            for raw_key, raw_value in raw_attrs.items():
                if alias in raw_key.lower() and raw_value and raw_value != "N/A":
                    value = raw_value
                    break
            if value:
                break
        harmonized[standard_field] = value

    return harmonized


def calculate_quality_score(record: dict, scoring_mode: str = "auto") -> dict:
    """Score datasets based on quality metrics. Returns Banked/Pending Review/Not Suitable."""
    scores = {}

    is_long_read = record.get("is_nanopore", False) or "PACBIO" in record.get("platform", "").upper()
    use_long_read_scoring = scoring_mode == "nanopore" or (scoring_mode == "auto" and is_long_read)

    if use_long_read_scoring:
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

    harmonized_count = record.get("harmonized_count", 0)
    if harmonized_count >= 8:
        scores["metadata_score"] = 25
    elif harmonized_count >= 5:
        scores["metadata_score"] = 20
    elif harmonized_count >= 3:
        scores["metadata_score"] = 15
    else:
        scores["metadata_score"] = 10

    harmonized = record.get("harmonized_fields", {})
    clinical_present = sum(1 for f in ["disease", "treatment", "age", "sex", "bmi"] if harmonized.get(f))
    if clinical_present >= 4:
        scores["clinical_score"] = 15
    elif clinical_present >= 2:
        scores["clinical_score"] = 10
    else:
        scores["clinical_score"] = 5

    scores["sample_score"] = 10 if record.get("is_fecal_sample") else 5
    scores["publication_score"] = 5 if record.get("pubmed_ids") else 0

    total_score = sum(scores.values())
    scores["total_score"] = total_score

    # Simplified 3-tier status
    if total_score >= 70:
        scores["vault_status"] = "Banked"
        scores["status_tier"] = 1
    elif total_score >= 55:
        scores["vault_status"] = "Pending Review"
        scores["status_tier"] = 2
    else:
        scores["vault_status"] = "Not Suitable"
        scores["status_tier"] = 3

    record.update(scores)
    return record


def create_results_dataframe(records: list[dict]) -> pd.DataFrame:
    """Create a pandas DataFrame from scored records."""
    if not records:
        return pd.DataFrame()

    df = pd.DataFrame(records)

    display_columns = [
        "run_accession", "run_url", "accession", "sra_url",
        "bioproject_id", "bioproject_url", "pubmed_ids",
        "title", "organism", "platform", "seq_type",
        "total_gb", "avg_read_length",
        "access_type", "vault_status", "disease_category",
        "harmonized_count", "is_fecal_sample", "gut_brain_relevant",
        "total_score"
    ]

    available_columns = [col for col in display_columns if col in df.columns]
    df = df[available_columns]

    if "total_score" in df.columns:
        df = df.sort_values("total_score", ascending=False)

    return df


def generate_one_line_summary(df: pd.DataFrame) -> str:
    """Generate a single sentence summary for quick understanding."""
    if df.empty:
        return "No datasets found. Try a different search."

    total = len(df)
    banked = len(df[df["vault_status"] == "Banked"]) if "vault_status" in df.columns else 0
    public = len(df[df["access_type"] == "Public"]) if "access_type" in df.columns else 0
    gut_brain = df["gut_brain_relevant"].sum() if "gut_brain_relevant" in df.columns else 0

    if banked > 0 and public > 0:
        return f"Found {total} datasets: {banked} are high-quality and {public} are publicly available for immediate use."
    elif banked > 0:
        return f"Found {total} datasets: {banked} meet quality thresholds for the vault."
    else:
        return f"Found {total} datasets, but none meet quality thresholds. Consider expanding your search."


def render_discovery_funnel(df: pd.DataFrame):
    """Render a visual funnel showing the discovery pipeline."""
    if df.empty:
        return

    total = len(df)
    banked = len(df[df["vault_status"] == "Banked"]) if "vault_status" in df.columns else 0
    pending = len(df[df["vault_status"] == "Pending Review"]) if "vault_status" in df.columns else 0
    public_banked = len(df[(df["vault_status"] == "Banked") & (df["access_type"] == "Public")]) if "vault_status" in df.columns and "access_type" in df.columns else 0

    fig = go.Figure(go.Funnel(
        y=["Discovered", "High Quality", "Public & Ready"],
        x=[total, banked, public_banked],
        textinfo="value+percent initial",
        marker={"color": ["#3498db", "#2ecc71", "#27ae60"]}
    ))

    fig.update_layout(
        title="Discovery Pipeline",
        height=250,
        margin=dict(l=20, r=20, t=40, b=20)
    )

    st.plotly_chart(fig, use_container_width=True)


def render_key_metrics(df: pd.DataFrame):
    """Render the 4 key metrics header."""
    if df.empty:
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Identified", 0)
        with col2:
            st.metric("Banked", 0)
        with col3:
            st.metric("Public", 0)
        with col4:
            st.metric("Gut-Brain", 0)
        return

    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric("Identified", len(df), help="Total datasets found")

    with col2:
        banked = len(df[df["vault_status"] == "Banked"]) if "vault_status" in df.columns else 0
        st.metric("Banked", banked, help="High-quality, ready for vault")

    with col3:
        public = len(df[df["access_type"] == "Public"]) if "access_type" in df.columns else 0
        st.metric("Public", public, help="Immediately downloadable")

    with col4:
        gut_brain = int(df["gut_brain_relevant"].sum()) if "gut_brain_relevant" in df.columns else 0
        st.metric("Gut-Brain", gut_brain, help="Relevant to gut-brain axis research")


def render_overview_tab(df: pd.DataFrame, records: list[dict]):
    """Render the Overview tab with summary, funnel, and insights."""
    # One-line summary
    summary = generate_one_line_summary(df)
    st.info(summary)

    if df.empty:
        return

    col1, col2 = st.columns([1, 1])

    with col1:
        render_discovery_funnel(df)

    with col2:
        # Quick insights
        st.markdown("### Quick Insights")

        total = len(df)
        banked = len(df[df["vault_status"] == "Banked"]) if "vault_status" in df.columns else 0
        public = len(df[df["access_type"] == "Public"]) if "access_type" in df.columns else 0
        gut_brain = int(df["gut_brain_relevant"].sum()) if "gut_brain_relevant" in df.columns else 0
        nanopore = len(df[df["seq_type"].str.contains("Nanopore", na=False)]) if "seq_type" in df.columns else 0

        insights = []
        if banked > 5:
            insights.append(f"**{banked} datasets ready** for bulk ingestion")
        elif banked > 0:
            insights.append(f"**{banked} datasets** meet quality thresholds")

        if gut_brain > 0:
            insights.append(f"**{gut_brain} datasets** relevant to gut-brain research")

        if nanopore > 0:
            insights.append(f"**{nanopore} long-read** datasets available")

        private = len(df[df["access_type"] == "Private (Requires Access)"]) if "access_type" in df.columns else 0
        if private > 0:
            insights.append(f"**{private} private datasets** could be pursued via partnership")

        if "disease_category" in df.columns:
            top = df["disease_category"].value_counts().head(1)
            if not top.empty and top.index[0] != "Unclassified":
                insights.append(f"Top category: **{top.index[0]}** ({top.values[0]} datasets)")

        for insight in insights:
            st.markdown(f"- {insight}")

        # Recommendation box
        st.markdown("### Recommendation")
        if banked > 5 and public > 3:
            st.success("Proceed with bulk ingestion - sufficient high-quality public data available.")
        elif banked > 0:
            st.warning("Review individual datasets before ingestion.")
        else:
            st.error("Expand search criteria - no datasets meet quality thresholds.")


def render_dataset_browser(df: pd.DataFrame, records: list[dict]):
    """Render the Dataset Browser tab with filterable table."""
    if df.empty:
        st.info("Run a search to see datasets.")
        return

    # Filters
    col1, col2, col3 = st.columns(3)

    with col1:
        if "vault_status" in df.columns:
            status_filter = st.multiselect(
                "Status",
                options=["Banked", "Pending Review", "Not Suitable"],
                default=["Banked", "Pending Review"]
            )
        else:
            status_filter = []

    with col2:
        if "access_type" in df.columns:
            access_filter = st.multiselect(
                "Access",
                options=df["access_type"].unique().tolist(),
                default=df["access_type"].unique().tolist()
            )
        else:
            access_filter = []

    with col3:
        fecal_only = st.checkbox("Fecal/gut samples only", value=False)

    # Apply filters
    filtered_df = df.copy()
    if "vault_status" in df.columns and status_filter:
        filtered_df = filtered_df[filtered_df["vault_status"].isin(status_filter)]
    if "access_type" in df.columns and access_filter:
        filtered_df = filtered_df[filtered_df["access_type"].isin(access_filter)]
    if fecal_only and "is_fecal_sample" in df.columns:
        filtered_df = filtered_df[filtered_df["is_fecal_sample"] == True]

    st.caption(f"Showing {len(filtered_df)} of {len(df)} datasets")

    # Display table
    st.dataframe(
        filtered_df,
        use_container_width=True,
        height=400,
        column_config={
            "run_accession": st.column_config.TextColumn("Run ID"),
            "run_url": st.column_config.LinkColumn("Link", display_text="View"),
            "vault_status": st.column_config.TextColumn("Status"),
            "access_type": st.column_config.TextColumn("Access"),
            "disease_category": st.column_config.TextColumn("Disease"),
            "total_gb": st.column_config.NumberColumn("Gb", format="%.2f"),
            "total_score": st.column_config.NumberColumn("Score"),
        }
    )

    # Export buttons
    col1, col2 = st.columns(2)
    with col1:
        st.download_button(
            "Download All Results (CSV)",
            filtered_df.to_csv(index=False),
            "discovery_results.csv",
            "text/csv"
        )
    with col2:
        banked_df = filtered_df[filtered_df["vault_status"] == "Banked"] if "vault_status" in filtered_df.columns else filtered_df
        if not banked_df.empty:
            st.download_button(
                "Download Banked Only (CSV)",
                banked_df.to_csv(index=False),
                "banked_datasets.csv",
                "text/csv"
            )

    # Verified source links
    st.markdown("### Verified NCBI Source Links")
    st.caption("Each dataset links to its official NCBI record. All accessions retrieved via E-utilities API.")

    for idx, row in filtered_df.head(10).iterrows():
        run = row.get("run_accession", "")
        run_url = row.get("run_url", "")
        bioproject = row.get("bioproject_id", "")
        bioproject_url = row.get("bioproject_url", "")
        pubmed = row.get("pubmed_ids", "")
        status = row.get("vault_status", "")
        title = str(row.get("title", ""))[:60]
        organism = row.get("organism", "N/A")
        disease = row.get("disease_category", "")

        links = []
        if run_url and run:
            links.append(f"[SRA: {run}]({run_url})")
        if bioproject_url and bioproject:
            links.append(f"[BioProject: {bioproject}]({bioproject_url})")
        if pubmed:
            for pmid in str(pubmed).split(",")[:2]:
                pmid = pmid.strip()
                if pmid:
                    links.append(f"[PubMed: {pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)")

        status_color = "green" if status == "Banked" else "orange" if status == "Pending Review" else "red"
        if links:
            st.markdown(f"**[{status}]** {title}...")
            st.markdown(f"*{organism}* | {disease} | {' | '.join(links)}")


def render_disease_cohorts(df: pd.DataFrame):
    """Render the Disease Cohorts tab."""
    if df.empty or "disease_category" not in df.columns:
        st.info("Run a search to see disease cohorts.")
        return

    col1, col2 = st.columns(2)

    with col1:
        # Disease distribution
        disease_counts = df["disease_category"].value_counts()
        fig = px.bar(
            x=disease_counts.values,
            y=disease_counts.index,
            orientation="h",
            title="Datasets by Disease Category",
            labels={"x": "Count", "y": "Category"},
            color=disease_counts.values,
            color_continuous_scale="Blues"
        )
        fig.update_layout(yaxis={"categoryorder": "total ascending"}, showlegend=False, height=350)
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        # Gut-brain relevance
        if "gut_brain_relevant" in df.columns:
            gb_count = int(df["gut_brain_relevant"].sum())
            other_count = len(df) - gb_count
            fig = px.pie(
                values=[gb_count, other_count],
                names=["Gut-Brain Relevant", "Other"],
                title="Gut-Brain Axis Relevance",
                color_discrete_sequence=["#9b59b6", "#bdc3c7"],
                hole=0.4
            )
            fig.update_layout(height=350)
            st.plotly_chart(fig, use_container_width=True)

    # Priority table
    st.markdown("### Cohort Prioritization")

    priority_data = []
    for category in df["disease_category"].unique():
        cat_df = df[df["disease_category"] == category]
        total = len(cat_df)
        banked = len(cat_df[cat_df["vault_status"] == "Banked"]) if "vault_status" in cat_df.columns else 0
        public = len(cat_df[cat_df["access_type"] == "Public"]) if "access_type" in cat_df.columns else 0

        if banked >= 3 and public >= 2:
            priority = "HIGH"
        elif banked >= 1:
            priority = "MEDIUM"
        else:
            priority = "LOW"

        priority_data.append({
            "Category": category,
            "Total": total,
            "Banked": banked,
            "Public": public,
            "Priority": priority
        })

    priority_df = pd.DataFrame(priority_data).sort_values("Banked", ascending=False)
    st.dataframe(priority_df, use_container_width=True, hide_index=True)


def render_data_quality(df: pd.DataFrame):
    """Render the Data Quality tab."""
    if df.empty:
        st.info("Run a search to see quality metrics.")
        return

    col1, col2 = st.columns(2)

    with col1:
        # Status distribution
        if "vault_status" in df.columns:
            status_counts = df["vault_status"].value_counts()
            fig = px.pie(
                values=status_counts.values,
                names=status_counts.index,
                title="Vault Status Distribution",
                color=status_counts.index,
                color_discrete_map={
                    "Banked": "#2ecc71",
                    "Pending Review": "#f39c12",
                    "Not Suitable": "#e74c3c"
                },
                hole=0.4
            )
            fig.update_layout(height=300)
            st.plotly_chart(fig, use_container_width=True)

    with col2:
        # Access type
        if "access_type" in df.columns:
            access_counts = df["access_type"].value_counts()
            fig = px.pie(
                values=access_counts.values,
                names=access_counts.index,
                title="Public vs Private",
                color=access_counts.index,
                color_discrete_map={
                    "Public": "#3498db",
                    "Private (Requires Access)": "#95a5a6"
                },
                hole=0.4
            )
            fig.update_layout(height=300)
            st.plotly_chart(fig, use_container_width=True)

    col1, col2 = st.columns(2)

    with col1:
        # Sequencing technology
        if "seq_type" in df.columns:
            tech_counts = df["seq_type"].value_counts()
            fig = px.pie(
                values=tech_counts.values,
                names=tech_counts.index,
                title="Sequencing Technology"
            )
            fig.update_layout(height=300)
            st.plotly_chart(fig, use_container_width=True)

    with col2:
        # Fecal samples
        if "is_fecal_sample" in df.columns:
            fecal_count = int(df["is_fecal_sample"].sum())
            other = len(df) - fecal_count
            fig = px.pie(
                values=[fecal_count, other],
                names=["Fecal/Gut", "Other"],
                title="Sample Type",
                color_discrete_sequence=["#27ae60", "#bdc3c7"]
            )
            fig.update_layout(height=300)
            st.plotly_chart(fig, use_container_width=True)

    # Metadata harmonization section
    st.markdown("### Metadata Harmonization")
    st.markdown(
        "Raw metadata from public repositories uses inconsistent labels. "
        "This tool maps disparate field names to a standardized ontology."
    )

    with st.expander("View Harmonization Ontology"):
        ontology_data = []
        for field, aliases in METADATA_ONTOLOGY.items():
            ontology_data.append({
                "Standard Field": field.title(),
                "Mapped Aliases": ", ".join(aliases[:4]) + ("..." if len(aliases) > 4 else "")
            })
        st.dataframe(pd.DataFrame(ontology_data), use_container_width=True, hide_index=True)

    # Scoring criteria
    st.markdown("### Quality Scoring Criteria")
    st.markdown("""
| Component | Max Points | What It Measures |
|-----------|------------|------------------|
| Sequencing Depth | 25 | Read count (short-read) or throughput in Gb (long-read) |
| Read Length | 25 | Longer reads enable better taxonomic resolution |
| Metadata Quality | 25 | Number of harmonized clinical/sample fields |
| Clinical Relevance | 15 | Presence of disease, treatment, demographics |
| Sample Type | 10 | Fecal/gut samples prioritized for microbiome research |
| Publication | 5 | Linked PubMed ID adds credibility |

**Status Thresholds:** Banked (70+) | Pending Review (55-69) | Not Suitable (<55)
    """)


def render_export_tab(df: pd.DataFrame, records: list[dict]):
    """Render the Export tab."""
    if df.empty:
        st.info("Run a search to export data.")
        return

    st.markdown("### Download Options")

    col1, col2, col3 = st.columns(3)

    with col1:
        st.download_button(
            "Full Results (CSV)",
            df.to_csv(index=False),
            "all_datasets.csv",
            "text/csv",
            help="All discovered datasets"
        )

    with col2:
        if "run_accession" in df.columns:
            banked_df = df[df["vault_status"] == "Banked"] if "vault_status" in df.columns else df
            accessions = banked_df["run_accession"].dropna().tolist()
            st.download_button(
                "Accession List (TXT)",
                "\n".join(accessions),
                "accessions.txt",
                "text/plain",
                help="SRR IDs for bulk download"
            )

    with col3:
        if records:
            banked_records = [r for r in records if r.get("vault_status") == "Banked"]
            if banked_records:
                export = {
                    "export_date": datetime.now().isoformat(),
                    "total": len(banked_records),
                    "datasets": [
                        {
                            "run": r.get("run_accession"),
                            "bioproject": r.get("bioproject_id"),
                            "organism": r.get("organism"),
                            "disease": r.get("disease_category"),
                            "gb": r.get("total_gb"),
                            "url": r.get("run_url")
                        }
                        for r in banked_records
                    ]
                }
                st.download_button(
                    "Banked JSON",
                    json.dumps(export, indent=2),
                    "banked.json",
                    "application/json",
                    help="Structured data for pipelines"
                )

    # Bulk download instructions
    st.markdown("### Bulk Download")
    st.code("""# Download using SRA Toolkit
prefetch --option-file accessions.txt
fasterq-dump --split-files *.sra""", language="bash")

    # Storage estimate
    if "total_gb" in df.columns:
        banked_df = df[df["vault_status"] == "Banked"] if "vault_status" in df.columns else df
        total_gb = banked_df["total_gb"].sum()

        st.markdown("### Storage Estimate")
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Banked Data Volume", f"{total_gb:.1f} Gb")
        with col2:
            st.metric("Est. S3 Cost", f"${total_gb * 0.023:.2f}/month")


def main():
    """Main application entry point."""
    st.set_page_config(
        page_title="Microbiome Discovery",
        page_icon=None,
        layout="wide"
    )

    st.title("Microbiome Dataset Discovery")
    st.caption("Automated discovery and curation for microbiome research")

    # Sidebar
    with st.sidebar:
        st.markdown("### Search")

        # Role-based presets
        role = st.selectbox("I am from...", list(ROLE_PRESETS.keys()))

        preset_options = list(ROLE_PRESETS[role].keys())
        preset = st.selectbox("Looking for...", ["Custom Query"] + preset_options)

        if preset == "Custom Query":
            query = st.text_area(
                "Query",
                value="fecal[All Fields] AND microbiome[All Fields]",
                height=80,
                help="NCBI Entrez syntax"
            )
        else:
            query = ROLE_PRESETS[role][preset]
            st.code(query, language=None)

        max_results = st.slider("Max results", 10, 100, 50, 10)

        st.markdown("---")

        scoring_mode = st.radio(
            "Scoring",
            ["auto", "nanopore", "illumina"],
            format_func=lambda x: {"auto": "Auto", "nanopore": "Long-Read", "illumina": "Short-Read"}[x],
            horizontal=True
        )

        search_button = st.button("Search", type="primary", use_container_width=True)

        st.markdown("---")
        st.caption(
            "Data from [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra)"
        )

    # Session state
    if "results_df" not in st.session_state:
        st.session_state.results_df = pd.DataFrame()
    if "records" not in st.session_state:
        st.session_state.records = []

    # Execute search
    if search_button:
        with st.spinner("Searching NCBI SRA..."):
            id_list = search_sra(query, max_results)

        if id_list:
            with st.spinner("Processing metadata..."):
                records = fetch_sra_metadata(id_list)
                scored_records = [calculate_quality_score(r, scoring_mode) for r in records]
                st.session_state.records = scored_records
                st.session_state.results_df = create_results_dataframe(scored_records)
        else:
            st.warning("No datasets found. Try adjusting your query.")

    df = st.session_state.results_df
    records = st.session_state.records

    # Key metrics
    render_key_metrics(df)

    # Tabs (simplified to 5)
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "Overview", "Dataset Browser", "Disease Cohorts", "Data Quality", "Export"
    ])

    with tab1:
        render_overview_tab(df, records)

    with tab2:
        render_dataset_browser(df, records)

    with tab3:
        render_disease_cohorts(df)

    with tab4:
        render_data_quality(df)

    with tab5:
        render_export_tab(df, records)

    # Footer with references and data sources
    st.markdown("---")
    st.markdown("### Data Sources & References")
    st.markdown("""
All data is retrieved in real-time from official NCBI repositories. Every accession number displayed is real and verifiable.

| Source | Description | Link |
|--------|-------------|------|
| NCBI SRA | Sequence Read Archive - primary source for all sequencing datasets | [ncbi.nlm.nih.gov/sra](https://www.ncbi.nlm.nih.gov/sra) |
| NCBI E-utilities | API used for programmatic data retrieval | [E-utilities Documentation](https://www.ncbi.nlm.nih.gov/books/NBK25501/) |
| NCBI BioProject | Study-level metadata and project information | [ncbi.nlm.nih.gov/bioproject](https://www.ncbi.nlm.nih.gov/bioproject/) |
| PubMed | Linked scientific publications | [pubmed.ncbi.nlm.nih.gov](https://pubmed.ncbi.nlm.nih.gov/) |

**Accession Types:**
- **SRR** (e.g., SRR12345678): Unique sequencing run - one per dataset
- **SRX** (e.g., SRX12345678): Experiment accession
- **PRJNA** (e.g., PRJNA123456): BioProject - groups related samples from a study
- **PMID** (e.g., 12345678): PubMed publication identifier
    """)

    st.caption(
        "Click any accession link to verify on the official NCBI website. "
        "Data quality scores are calculated based on sequencing depth, read length, "
        "metadata completeness, and clinical relevance."
    )


if __name__ == "__main__":
    main()
