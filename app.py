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

# Disease categories for prioritization (gut-brain axis focus)
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

# Study type classification
STUDY_TYPES = {
    "16S Amplicon": ["16s", "amplicon", "v3-v4", "v4", "rrna", "16s rrna"],
    "Shotgun Metagenomics": ["wgs", "shotgun", "metagenome", "metagenomic", "whole genome"],
    "Isolate Genome": ["isolate", "pure culture", "single strain", "genome assembly"],
    "Metatranscriptomics": ["rna-seq", "metatranscript", "transcriptome"],
    "Clinical Trial": ["randomized", "placebo", "intervention", "clinical trial", "rct"]
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

    # Classify disease relevance from title/abstract
    searchable_text = " ".join([
        record.get("title", ""),
        record.get("study_title", ""),
        record.get("study_abstract", "")
    ]).lower()

    # Disease classification
    detected_diseases = []
    for category, keywords in DISEASE_CATEGORIES.items():
        if any(kw in searchable_text for kw in keywords):
            detected_diseases.append(category)
    record["disease_categories"] = detected_diseases
    record["disease_category"] = detected_diseases[0] if detected_diseases else "Unclassified"
    record["gut_brain_relevant"] = "Gut-Brain Axis" in detected_diseases or "Pain Conditions" in detected_diseases

    # Study type classification
    detected_study_types = []
    library_strategy = record.get("library_strategy", "").lower()
    for study_type, keywords in STUDY_TYPES.items():
        if any(kw in searchable_text or kw in library_strategy for kw in keywords):
            detected_study_types.append(study_type)
    record["study_types"] = detected_study_types
    record["primary_study_type"] = detected_study_types[0] if detected_study_types else "Other"

    # Estimate download size (rough: 1 read ~ 500 bytes compressed)
    record["estimated_size_gb"] = (record.get("total_bases", 0) * 0.3) / 1e9  # ~30% compression

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


def render_executive_summary(df: pd.DataFrame):
    """Render executive summary for stakeholders."""
    if df.empty:
        return

    st.markdown("## Executive Summary")

    # Key insight box
    total = len(df)
    banked = len(df[df["ingestion_decision"] == "AUTO-INGEST"]) if "ingestion_decision" in df.columns else 0
    public = len(df[df["access_type"] == "Public"]) if "access_type" in df.columns else 0
    nanopore = len(df[df["seq_type"].str.contains("Nanopore", na=False)]) if "seq_type" in df.columns else 0

    col1, col2 = st.columns([2, 1])

    with col1:
        st.success(f"""
        **Discovery Complete** - {total} datasets identified

        - **{banked} datasets** ({banked/total*100:.0f}%) meet quality thresholds for immediate ingestion
        - **{public} datasets** are publicly accessible (no restrictions)
        - **{nanopore} datasets** use Oxford Nanopore long-read sequencing
        """)

    with col2:
        # Progress indicator
        progress = banked / total if total > 0 else 0
        st.metric("Pipeline Ready", f"{progress*100:.0f}%")
        st.progress(progress)

    st.markdown("---")


def render_summary_header(df: pd.DataFrame):
    """Render the summary header with live statistics."""
    st.markdown("## Live Discovery Status")
    st.caption(
        "Data Source: [NCBI Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) | "
        "[E-utilities API Documentation](https://www.ncbi.nlm.nih.gov/books/NBK25501/)"
    )

    if df.empty:
        col1, col2, col3, col4, col5, col6 = st.columns(6)
        with col1:
            st.metric("Identified", 0)
        with col2:
            st.metric("Banked", 0)
        with col3:
            st.metric("Pending Review", 0)
        with col4:
            st.metric("Public", 0)
        with col5:
            st.metric("Restricted", 0)
        with col6:
            st.metric("Total Gb", "0.0")
        return

    col1, col2, col3, col4, col5, col6 = st.columns(6)

    with col1:
        st.metric("Identified", len(df), help="Total datasets found matching query")

    with col2:
        if "ingestion_decision" in df.columns:
            banked = len(df[df["ingestion_decision"] == "AUTO-INGEST"])
            st.metric("Banked", banked,
                     delta=f"{banked/len(df)*100:.0f}%" if len(df) > 0 else "0%",
                     help="High-quality datasets ready for ingestion")

    with col3:
        if "ingestion_decision" in df.columns:
            review = len(df[df["ingestion_decision"] == "REVIEW"])
            st.metric("Pending Review", review, help="Datasets requiring manual review")

    with col4:
        if "access_type" in df.columns:
            public = len(df[df["access_type"] == "Public"])
            st.metric("Public", public, help="Publicly accessible datasets")

    with col5:
        if "access_type" in df.columns:
            restricted = len(df[df["access_type"].str.contains("Restricted", na=False)])
            st.metric("Restricted", restricted, help="Controlled access (dbGaP, etc.)")

    with col6:
        if "total_gb" in df.columns:
            total_gb = df["total_gb"].sum()
            st.metric("Total Gb", f"{total_gb:.1f}", help="Total sequencing data volume")


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
            st.plotly_chart(fig, width="stretch")

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
            st.plotly_chart(fig, width="stretch")

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
            st.plotly_chart(fig, width="stretch")

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
            st.plotly_chart(fig, width="stretch")

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
            st.plotly_chart(fig, width="stretch")

    with col2:
        if "is_fecal_sample" in df.columns:
            fecal_counts = df["is_fecal_sample"].value_counts()
            fig = px.pie(
                values=fecal_counts.values,
                names=["Fecal/Gut" if x else "Other" for x in fecal_counts.index],
                title="Sample Type (Gut Relevance)"
            )
            st.plotly_chart(fig, width="stretch")


def render_organism_breakdown(df: pd.DataFrame):
    """Show organism/taxonomy breakdown for microbiologists."""
    st.subheader("Organism Distribution")
    st.markdown("Taxonomic breakdown of organisms in the discovered datasets.")

    if df.empty or "organism" not in df.columns:
        st.info("No organism data available.")
        return

    col1, col2 = st.columns(2)

    with col1:
        # Top organisms
        organism_counts = df["organism"].value_counts().head(10)
        if not organism_counts.empty:
            fig = px.bar(
                x=organism_counts.values,
                y=organism_counts.index,
                orientation="h",
                title="Top 10 Organisms",
                labels={"x": "Dataset Count", "y": "Organism"}
            )
            fig.update_layout(yaxis={"categoryorder": "total ascending"})
            st.plotly_chart(fig, width="stretch")

    with col2:
        # Organism categories
        def categorize_organism(org):
            org_lower = str(org).lower()
            if "human" in org_lower or "homo" in org_lower:
                return "Human microbiome"
            elif "mouse" in org_lower or "mus" in org_lower:
                return "Mouse microbiome"
            elif "gut" in org_lower or "fecal" in org_lower or "intestin" in org_lower:
                return "Gut-associated"
            elif "metagenome" in org_lower:
                return "Metagenome"
            elif "bacteria" in org_lower:
                return "Bacterial isolate"
            else:
                return "Other"

        df_copy = df.copy()
        df_copy["organism_category"] = df_copy["organism"].apply(categorize_organism)
        cat_counts = df_copy["organism_category"].value_counts()

        fig = px.pie(
            values=cat_counts.values,
            names=cat_counts.index,
            title="Sample Categories"
        )
        st.plotly_chart(fig, width="stretch")


def render_sample_type_analysis(df: pd.DataFrame, records: list[dict]):
    """Show sample type breakdown."""
    st.subheader("Sample Type Analysis")
    st.markdown("Distribution of sample types relevant for microbiome research.")

    if not records:
        st.info("No sample data available.")
        return

    # Extract sample types from harmonized fields
    sample_types = []
    body_sites = []

    for rec in records:
        harmonized = rec.get("harmonized_fields", {})
        sample_type = harmonized.get("sample_type", "")
        if sample_type:
            sample_types.append(sample_type.lower())

        # Check raw attributes for body site info
        attrs = rec.get("sample_attributes", {})
        for key, value in attrs.items():
            if any(term in key.lower() for term in ["body_site", "tissue", "isolation_source"]):
                if value and value != "N/A":
                    body_sites.append(value.lower())

    col1, col2 = st.columns(2)

    with col1:
        if sample_types:
            type_counts = pd.Series(sample_types).value_counts().head(10)
            fig = px.bar(
                x=type_counts.values,
                y=type_counts.index,
                orientation="h",
                title="Sample Types",
                labels={"x": "Count", "y": "Type"}
            )
            fig.update_layout(yaxis={"categoryorder": "total ascending"})
            st.plotly_chart(fig, width="stretch")
        else:
            st.info("No sample type metadata available.")

    with col2:
        # Fecal vs other breakdown
        fecal_count = len(df[df["is_fecal_sample"] == True]) if "is_fecal_sample" in df.columns else 0
        other_count = len(df) - fecal_count

        fig = px.pie(
            values=[fecal_count, other_count],
            names=["Fecal/Gut", "Other"],
            title="Gut Microbiome Relevance",
            color_discrete_sequence=["#2ecc71", "#95a5a6"]
        )
        st.plotly_chart(fig, width="stretch")


def render_disease_prioritization(df: pd.DataFrame, records: list[dict]):
    """Show disease prioritization for strategic planning."""
    st.subheader("Disease Cohort Prioritization")
    st.markdown(
        "Strategic overview of disease categories represented in discovered datasets. "
        "Helps prioritize which cohorts to pursue for data acquisition."
    )

    if df.empty:
        st.info("No data available.")
        return

    col1, col2 = st.columns(2)

    with col1:
        # Disease category distribution
        if "disease_category" in df.columns:
            disease_counts = df["disease_category"].value_counts()
            fig = px.bar(
                x=disease_counts.values,
                y=disease_counts.index,
                orientation="h",
                title="Datasets by Disease Category",
                labels={"x": "Dataset Count", "y": "Disease Category"},
                color=disease_counts.values,
                color_continuous_scale="Blues"
            )
            fig.update_layout(yaxis={"categoryorder": "total ascending"}, showlegend=False)
            st.plotly_chart(fig, width="stretch")

    with col2:
        # Gut-brain axis relevance
        if "gut_brain_relevant" in df.columns:
            gb_count = df["gut_brain_relevant"].sum()
            other_count = len(df) - gb_count
            fig = px.pie(
                values=[gb_count, other_count],
                names=["Gut-Brain Relevant", "Other"],
                title="Gut-Brain Axis Relevance",
                color_discrete_sequence=["#9b59b6", "#95a5a6"]
            )
            st.plotly_chart(fig, width="stretch")

    # Priority recommendations
    st.markdown("### Prioritization Recommendations")

    if "disease_category" in df.columns and "ingestion_decision" in df.columns:
        priority_data = []
        for category in df["disease_category"].unique():
            cat_df = df[df["disease_category"] == category]
            total = len(cat_df)
            banked = len(cat_df[cat_df["ingestion_decision"] == "AUTO-INGEST"])
            public = len(cat_df[cat_df["access_type"] == "Public"]) if "access_type" in cat_df.columns else 0
            total_gb = cat_df["total_gb"].sum() if "total_gb" in cat_df.columns else 0

            priority_data.append({
                "Disease Category": category,
                "Total Datasets": total,
                "High Quality": banked,
                "Public Access": public,
                "Data Volume (Gb)": f"{total_gb:.1f}",
                "Recommendation": "HIGH PRIORITY" if banked >= 3 and public >= 2 else "MEDIUM" if banked >= 1 else "LOW"
            })

        priority_df = pd.DataFrame(priority_data).sort_values("High Quality", ascending=False)
        st.dataframe(priority_df, width="stretch")


def render_study_type_analysis(df: pd.DataFrame):
    """Show study type classification."""
    st.subheader("Study Type Classification")
    st.markdown("Breakdown by sequencing approach (16S amplicon vs shotgun metagenomics vs isolate genomes).")

    if df.empty or "primary_study_type" not in df.columns:
        st.info("No study type data available.")
        return

    col1, col2 = st.columns(2)

    with col1:
        type_counts = df["primary_study_type"].value_counts()
        fig = px.pie(
            values=type_counts.values,
            names=type_counts.index,
            title="Study Types"
        )
        st.plotly_chart(fig, width="stretch")

    with col2:
        # Cross-tab: study type vs quality
        if "quality_grade" in df.columns:
            cross_tab = pd.crosstab(df["primary_study_type"], df["quality_grade"])
            fig = px.bar(
                cross_tab,
                barmode="stack",
                title="Study Type by Quality Grade",
                color_discrete_map={"A": "#2ecc71", "B": "#3498db", "C": "#f39c12", "D": "#e74c3c"}
            )
            st.plotly_chart(fig, width="stretch")


def render_actionable_insights(df: pd.DataFrame, records: list[dict]):
    """Generate actionable insights for stakeholders."""
    st.subheader("Actionable Insights")
    st.markdown("Auto-generated recommendations based on discovery results.")

    if df.empty:
        st.info("Run a search to generate insights.")
        return

    insights = []

    # Total stats
    total = len(df)
    banked = len(df[df["ingestion_decision"] == "AUTO-INGEST"]) if "ingestion_decision" in df.columns else 0
    public = len(df[df["access_type"] == "Public"]) if "access_type" in df.columns else 0

    insights.append(f"Discovered **{total} datasets** matching search criteria.")
    insights.append(f"**{banked} datasets** ({banked/total*100:.0f}%) meet quality thresholds for immediate ingestion.")
    insights.append(f"**{public} datasets** are publicly accessible without restrictions.")

    # Gut-brain specific
    if "gut_brain_relevant" in df.columns:
        gb_count = df["gut_brain_relevant"].sum()
        if gb_count > 0:
            insights.append(f"**{gb_count} datasets** are relevant to gut-brain axis research.")

    # Top disease recommendation
    if "disease_category" in df.columns:
        top_disease = df["disease_category"].value_counts().head(1)
        if not top_disease.empty:
            insights.append(f"Most represented category: **{top_disease.index[0]}** ({top_disease.values[0]} datasets).")

    # Data volume
    if "total_gb" in df.columns:
        total_gb = df["total_gb"].sum()
        insights.append(f"Total data volume: **{total_gb:.1f} Gb** across all datasets.")

    # Nanopore availability
    if "seq_type" in df.columns:
        nanopore = len(df[df["seq_type"].str.contains("Nanopore", na=False)])
        if nanopore > 0:
            insights.append(f"**{nanopore} datasets** use Oxford Nanopore long-read sequencing.")

    # Display insights
    for insight in insights:
        st.markdown(f"- {insight}")

    # Specific recommendations
    st.markdown("### Recommendations")

    recommendations = []

    if banked > 5:
        recommendations.append("**Proceed with bulk ingestion** - sufficient high-quality datasets available.")
    elif banked > 0:
        recommendations.append("**Selective ingestion recommended** - review individual datasets before ingestion.")
    else:
        recommendations.append("**Expand search criteria** - current results do not meet quality thresholds.")

    if "gut_brain_relevant" in df.columns and df["gut_brain_relevant"].sum() > 3:
        recommendations.append("**Gut-brain cohort opportunity** - multiple relevant datasets identified.")

    restricted = len(df[df["access_type"].str.contains("Restricted", na=False)]) if "access_type" in df.columns else 0
    if restricted > 0:
        recommendations.append(f"**{restricted} restricted datasets** require dbGaP/controlled access applications.")

    for rec in recommendations:
        st.success(rec)


def render_ml_export(df: pd.DataFrame, records: list[dict]):
    """Export data for ML pipelines."""
    st.subheader("ML Pipeline Export")
    st.markdown("Export data in formats ready for downstream ML pipelines and bulk download scripts.")

    if df.empty:
        st.info("No data to export.")
        return

    col1, col2, col3 = st.columns(3)

    with col1:
        # CSV export
        csv_data = df.to_csv(index=False)
        st.download_button(
            "Download CSV",
            csv_data,
            "microbiome_datasets.csv",
            "text/csv",
            help="Full dataset with all metadata"
        )

    with col2:
        # Accession list for bulk download
        if "run_accession" in df.columns:
            accessions = df["run_accession"].dropna().tolist()
            accession_text = "\n".join(accessions)
            st.download_button(
                "Download Accession List",
                accession_text,
                "accessions.txt",
                "text/plain",
                help="SRR accessions for sra-toolkit prefetch/fasterq-dump"
            )

    with col3:
        # JSON export for ML
        if records:
            # Filter to banked only
            banked_records = [r for r in records if r.get("ingestion_decision") == "AUTO-INGEST"]
            if banked_records:
                json_export = {
                    "export_date": datetime.now().isoformat(),
                    "total_datasets": len(banked_records),
                    "datasets": [
                        {
                            "run_accession": r.get("run_accession"),
                            "experiment_accession": r.get("accession"),
                            "bioproject": r.get("bioproject_id"),
                            "organism": r.get("organism"),
                            "platform": r.get("platform"),
                            "total_bases": r.get("total_bases"),
                            "total_gb": r.get("total_gb"),
                            "quality_score": r.get("total_score"),
                            "disease_category": r.get("disease_category"),
                            "study_type": r.get("primary_study_type"),
                            "harmonized_metadata": r.get("harmonized_fields", {}),
                            "sra_url": r.get("sra_url"),
                            "run_url": r.get("run_url")
                        }
                        for r in banked_records
                    ]
                }
                st.download_button(
                    "Download JSON (Banked)",
                    json.dumps(json_export, indent=2),
                    "ingestion_manifest.json",
                    "application/json",
                    help="Structured JSON for ML pipeline ingestion"
                )

    # Infrastructure estimates
    st.markdown("### Infrastructure Estimates")

    if "estimated_size_gb" in df.columns:
        total_size = df["estimated_size_gb"].sum()
        banked_df = df[df["ingestion_decision"] == "AUTO-INGEST"] if "ingestion_decision" in df.columns else df
        banked_size = banked_df["estimated_size_gb"].sum() if "estimated_size_gb" in banked_df.columns else 0

        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Total Est. Download", f"{total_size:.1f} Gb")
        with col2:
            st.metric("Banked Est. Download", f"{banked_size:.1f} Gb")
        with col3:
            # S3 cost estimate (~$0.023/GB/month for standard)
            monthly_cost = banked_size * 0.023
            st.metric("Est. S3 Storage Cost", f"${monthly_cost:.2f}/mo")

    # Bulk download command
    st.markdown("### Bulk Download Command")
    st.code("""# Using SRA Toolkit
prefetch --option-file accessions.txt
fasterq-dump --split-files *.sra

# Or using parallel download
cat accessions.txt | xargs -P 4 -I {} prefetch {}""", language="bash")


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
        st.dataframe(pd.DataFrame(ontology_df), width="stretch")

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

        st.dataframe(coverage_df, width="stretch")


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
        st.markdown("**Quick Search**")
        preset = st.selectbox(
            "Preset Queries",
            [
                "Custom",
                "Gut-Brain Axis Studies",
                "Depression & Microbiome",
                "Anxiety & Gut Microbiome",
                "IBS / IBD Studies",
                "Pain & Microbiome",
                "Nanopore Fecal Studies",
                "Human Gut Microbiome",
                "Clinical Stool Studies",
                "Probiotic Trials",
                "Prebiotic Studies",
                "16S Amplicon (Any Platform)",
                "Shotgun Metagenomics"
            ]
        )

        preset_queries = {
            "Custom": "",
            "Gut-Brain Axis Studies": "(gut-brain[All Fields] OR gut brain axis[All Fields]) AND microbiome[All Fields]",
            "Depression & Microbiome": "(depression[All Fields] OR depressive[All Fields]) AND (gut[All Fields] OR fecal[All Fields]) AND microbiome[All Fields]",
            "Anxiety & Gut Microbiome": "anxiety[All Fields] AND (gut[All Fields] OR fecal[All Fields]) AND microbiome[All Fields]",
            "IBS / IBD Studies": "(IBS[All Fields] OR IBD[All Fields] OR \"irritable bowel\"[All Fields]) AND microbiome[All Fields]",
            "Pain & Microbiome": "(pain[All Fields] OR nociception[All Fields]) AND (gut[All Fields] OR microbiome[All Fields])",
            "Nanopore Fecal Studies": "fecal[All Fields] AND Oxford Nanopore[Platform]",
            "Human Gut Microbiome": "gut microbiome[All Fields] AND human[Organism]",
            "Clinical Stool Studies": "stool[All Fields] AND clinical[All Fields] AND microbiome[All Fields]",
            "Probiotic Trials": "probiotic[All Fields] AND (clinical trial[All Fields] OR randomized[All Fields]) AND gut[All Fields]",
            "Prebiotic Studies": "prebiotic[All Fields] AND (gut[All Fields] OR microbiome[All Fields])",
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

        search_button = st.button("Launch Discovery Agent", type="primary", width="stretch")

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

    # Show summary header (live metrics)
    render_summary_header(df)

    if not df.empty:
        # Executive summary for stakeholders
        render_executive_summary(df)

        # Tabs
        tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9, tab10 = st.tabs([
            "Dataset Queue", "Quality Analytics", "Disease Priority", "Study Types",
            "Insights", "ML Export", "Organisms", "Sample Types", "Metadata", "Nanopore"
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
                width="stretch",
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

            # Verified source links with NCBI references
            st.markdown("### Verified NCBI Source Links")
            st.markdown(
                "Each dataset below links to its official NCBI record. "
                "All accessions are verified via [NCBI E-utilities API](https://www.ncbi.nlm.nih.gov/books/NBK25501/)."
            )

            for idx, row in filtered_df.head(10).iterrows():
                run = row.get("run_accession", "")
                run_url = row.get("run_url", "")
                bioproject = row.get("bioproject_id", "")
                bioproject_url = row.get("bioproject_url", "")
                pubmed = row.get("pubmed_ids", "")
                decision = row.get("ingestion_decision", "")
                title = str(row.get("title", ""))[:50]
                organism = row.get("organism", "")

                links = []
                if run_url and run:
                    links.append(f"[SRA Run: {run}]({run_url})")
                if bioproject_url and bioproject:
                    links.append(f"[BioProject: {bioproject}]({bioproject_url})")
                if pubmed:
                    for pmid in str(pubmed).split(","):
                        pmid = pmid.strip()
                        if pmid:
                            links.append(f"[PubMed: {pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)")

                status_icon = "BANKED" if decision == "AUTO-INGEST" else "REVIEW" if decision == "REVIEW" else "REJECT"
                if links:
                    st.markdown(f"**[{status_icon}]** {title}... (*{organism}*)")
                    st.markdown(f"  {' | '.join(links)}")

        with tab2:
            st.subheader("Quality Analytics")
            render_quality_charts(df)

        with tab3:
            render_disease_prioritization(df, records)

        with tab4:
            render_study_type_analysis(df)

        with tab5:
            render_actionable_insights(df, records)

        with tab6:
            render_ml_export(df, records)

        with tab7:
            render_organism_breakdown(df)

        with tab8:
            render_sample_type_analysis(df, records)

        with tab9:
            render_metadata_harmonization(df, records)

        with tab10:
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
                    st.plotly_chart(fig, width="stretch")

                st.dataframe(long_read_df, width="stretch", height=300)
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

    # Footer with data source attribution
    st.markdown("---")
    st.markdown(
        "**Data Sources & References**\n\n"
        "- [NCBI Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) - Primary data source\n"
        "- [NCBI E-utilities API](https://www.ncbi.nlm.nih.gov/books/NBK25501/) - Data retrieval\n"
        "- [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/) - Study metadata\n"
        "- [PubMed](https://pubmed.ncbi.nlm.nih.gov/) - Linked publications\n"
    )
    st.caption(
        "All dataset accessions (SRR, SRX, PRJNA) are real identifiers retrieved from NCBI. "
        "Click any link to verify on the official NCBI website."
    )


if __name__ == "__main__":
    main()
