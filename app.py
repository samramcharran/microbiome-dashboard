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
import re


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

# Disease burden weights for Mission Control prioritization
DISEASE_BURDEN_WEIGHTS = {
    "Gut-Brain Axis": 10,
    "Pain Conditions": 9,
    "GI Disorders": 8,
    "Metabolic": 7,
    "Immune/Inflammatory": 6,
    "Infectious": 5,
    "Healthy/Control": 3,
    "Unclassified": 1
}

# Leadership dashboard views (not searches - these control what metrics to show)
LEADERSHIP_VIEWS = [
    "Banking Progress",
    "Disease Prioritization",
    "Data Quality",
    "Access Status",
    "Infrastructure",
]

# Researcher search options - organized by category
RESEARCHER_SEARCHES = {
    # By Condition
    "Mental Health": "(depression[All Fields] OR anxiety[All Fields] OR stress[All Fields] OR autism[All Fields] OR parkinson[All Fields]) AND gut[All Fields] AND microbiome[All Fields]",
    "Pain Conditions": "(pain[All Fields] OR fibromyalgia[All Fields] OR chronic pain[All Fields] OR neuropathy[All Fields]) AND microbiome[All Fields]",
    "Digestive Health": "(IBS[All Fields] OR IBD[All Fields] OR crohn[All Fields] OR colitis[All Fields]) AND microbiome[All Fields]",
    "Metabolic Health": "(obesity[All Fields] OR diabetes[All Fields] OR metabolic[All Fields]) AND gut[All Fields] AND microbiome[All Fields]",
    # By Data Type
    "Long-Read Sequencing": "(fecal[All Fields] OR stool[All Fields]) AND (Oxford Nanopore[Platform] OR PacBio[Platform])",
    "Shotgun Metagenomics": "(shotgun[All Fields] OR metagenome[All Fields] OR WGS[All Fields]) AND fecal[All Fields] AND human[Organism]",
    "16S Amplicon": "16S[All Fields] AND (fecal[All Fields] OR stool[All Fields]) AND human[Organism]",
    # By Quality & Access
    "High-Quality Only": "fecal[All Fields] AND microbiome[All Fields] AND human[Organism]",  # Filtered post-search
    "Public Datasets": "fecal[All Fields] AND microbiome[All Fields] AND human[Organism]",  # Filtered post-search
    # By Study Design
    "Large Cohorts": "stool[All Fields] AND cohort[All Fields] AND microbiome[All Fields]",
    "Clinical Trials": "(clinical trial[All Fields] OR randomized[All Fields] OR intervention[All Fields]) AND gut[All Fields] AND microbiome[All Fields]",
    # Browse All
    "All Available Data": "fecal[All Fields] AND microbiome[All Fields] AND human[Organism]",
}

# Roles
ROLES = ["Researcher", "Leadership"]

# Microbiome loading messages
LOADING_MESSAGES = [
    "Exploring the gut microbiome...",
    "Analyzing microbial communities...",
    "Discovering beneficial bacteria...",
    "Mapping the microbiome landscape...",
    "Searching for microbial diversity...",
]


def render_loading_animation(message: str = "Loading..."):
    """Render a microbiome-themed loading animation."""
    # SVG animation of bacteria/microbes
    # Note: CSS curly braces are escaped as {{ }} for Python format string
    microbe_svg = """
    <div style="text-align: center; padding: 40px;">
        <svg width="120" height="120" viewBox="0 0 120 120">
            <style>
                .microbe {{ animation: float 2s ease-in-out infinite; }}
                .microbe1 {{ animation-delay: 0s; }}
                .microbe2 {{ animation-delay: 0.3s; }}
                .microbe3 {{ animation-delay: 0.6s; }}
                .microbe4 {{ animation-delay: 0.9s; }}
                .microbe5 {{ animation-delay: 1.2s; }}
                @keyframes float {{
                    0%, 100% {{ transform: translateY(0px); opacity: 0.7; }}
                    50% {{ transform: translateY(-10px); opacity: 1; }}
                }}
            </style>
            <!-- Bacteria shapes -->
            <ellipse class="microbe microbe1" cx="30" cy="60" rx="12" ry="8" fill="#4CAF50"/>
            <ellipse class="microbe microbe2" cx="60" cy="40" rx="10" ry="6" fill="#8BC34A"/>
            <ellipse class="microbe microbe3" cx="90" cy="55" rx="11" ry="7" fill="#66BB6A"/>
            <ellipse class="microbe microbe4" cx="45" cy="85" rx="9" ry="6" fill="#81C784"/>
            <ellipse class="microbe microbe5" cx="75" cy="80" rx="10" ry="7" fill="#A5D6A7"/>
            <!-- Small dots for flagella effect -->
            <circle class="microbe microbe1" cx="18" cy="60" r="3" fill="#4CAF50"/>
            <circle class="microbe microbe2" cx="50" cy="38" r="2" fill="#8BC34A"/>
            <circle class="microbe microbe3" cx="101" cy="53" r="3" fill="#66BB6A"/>
        </svg>
        <p style="color: #666; font-size: 16px; margin-top: 20px;">{message}</p>
    </div>
    """
    return microbe_svg.format(message=message)


def parse_flexible_date(date_str: str) -> Optional[datetime]:
    """Parse dates in multiple formats: YYYY-MM-DD, YYYY-MM, YYYY."""
    if not date_str or date_str == "N/A":
        return None

    date_str = str(date_str).strip()

    # Try various date formats
    formats = [
        "%Y-%m-%d",
        "%Y-%m",
        "%Y",
        "%d-%m-%Y",
        "%m-%d-%Y",
        "%Y/%m/%d",
        "%d/%m/%Y",
        "%m/%d/%Y",
        "%b %Y",
        "%B %Y",
        "%Y-%m-%dT%H:%M:%S",
        "%Y-%m-%dT%H:%M:%SZ",
    ]

    for fmt in formats:
        try:
            return datetime.strptime(date_str, fmt)
        except ValueError:
            continue

    # Try to extract year if nothing else works
    year_match = re.search(r'\b(19|20)\d{2}\b', date_str)
    if year_match:
        try:
            return datetime.strptime(year_match.group(), "%Y")
        except ValueError:
            pass

    return None


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


def render_temporal_trend(df: pd.DataFrame, records: list[dict]):
    """Render a bar chart showing datasets over time."""
    st.markdown("### Temporal Trend")

    if not records:
        st.caption("Run a search to see temporal distribution of datasets.")
        return

    dates = []
    for record in records:
        harmonized = record.get("harmonized_fields", {})
        date_str = harmonized.get("collection_date")
        if date_str:
            parsed_date = parse_flexible_date(date_str)
            if parsed_date:
                dates.append(parsed_date)

    if not dates:
        st.caption("No collection date metadata available for these datasets.")
        return

    # Aggregate by year
    years = [d.year for d in dates]
    year_counts = pd.Series(years).value_counts().sort_index()

    fig = px.bar(
        x=year_counts.index,
        y=year_counts.values,
        labels={"x": "Year", "y": "Datasets"},
        title=f"Datasets by Collection Year ({len(dates)} with dates)",
        color_discrete_sequence=["#3498db"]
    )
    fig.update_layout(height=300, margin=dict(l=20, r=20, t=40, b=20))
    st.plotly_chart(fig, use_container_width=True)


def render_platform_breakdown(df: pd.DataFrame):
    """Render a breakdown of sequencing platforms and data volume."""
    st.markdown("### Platform & Data Volume")

    if df.empty:
        st.caption("Run a search to see platform breakdown.")
        return

    col1, col2 = st.columns(2)

    with col1:
        # Platform pie chart
        if "seq_type" in df.columns:
            platform_counts = df["seq_type"].value_counts()
            fig = px.pie(
                values=platform_counts.values,
                names=platform_counts.index,
                title="Sequencing Technology",
                color_discrete_sequence=px.colors.qualitative.Set2
            )
            fig.update_layout(height=280, margin=dict(l=10, r=10, t=40, b=10))
            st.plotly_chart(fig, use_container_width=True)

    with col2:
        # Data volume by status
        if "vault_status" in df.columns and "total_gb" in df.columns:
            volume_by_status = df.groupby("vault_status")["total_gb"].sum().reset_index()
            fig = px.bar(
                volume_by_status,
                x="vault_status",
                y="total_gb",
                title="Data Volume by Status (Gb)",
                color="vault_status",
                color_discrete_map={
                    "Banked": "#2ecc71",
                    "Pending Review": "#f39c12",
                    "Not Suitable": "#e74c3c"
                }
            )
            fig.update_layout(height=280, margin=dict(l=10, r=10, t=40, b=10), showlegend=False)
            st.plotly_chart(fig, use_container_width=True)


def render_comparison_view(selected_accessions: list[str], records: list[dict]):
    """Render side-by-side comparison of selected datasets."""
    if not selected_accessions:
        st.info("Select datasets to compare from the browser above.")
        return

    # Find matching records
    selected_records = []
    for acc in selected_accessions:
        for record in records:
            if record.get("run_accession") == acc or record.get("accession") == acc:
                selected_records.append(record)
                break

    if not selected_records:
        st.warning("Could not find selected datasets.")
        return

    # Comparison fields
    comparison_fields = [
        ("Accession", "accession"),
        ("Run ID", "run_accession"),
        ("Organism", "organism"),
        ("Platform", "platform"),
        ("Seq Type", "seq_type"),
        ("Total Gb", "total_gb"),
        ("Quality Score", "total_score"),
        ("Status", "vault_status"),
        ("Disease Category", "disease_category"),
        ("Access Type", "access_type"),
        ("Fecal Sample", "is_fecal_sample"),
        ("Gut-Brain Relevant", "gut_brain_relevant"),
    ]

    # Create columns for side-by-side view
    cols = st.columns(len(selected_records))

    for idx, (col, record) in enumerate(zip(cols, selected_records)):
        with col:
            st.markdown(f"**Dataset {idx + 1}**")
            for label, field in comparison_fields:
                value = record.get(field, "N/A")
                if field == "total_gb" and isinstance(value, (int, float)):
                    value = f"{value:.2f}"
                elif isinstance(value, bool):
                    value = "Yes" if value else "No"
                st.markdown(f"**{label}:** {value}")

            # Links
            run_url = record.get("run_url", "")
            if run_url:
                st.markdown(f"[View on NCBI]({run_url})")


def render_disease_burden_chart(df: pd.DataFrame):
    """Render disease burden priority matrix visualization."""
    if df.empty or "disease_category" not in df.columns:
        return

    # Calculate weighted scores
    disease_data = []
    for category in df["disease_category"].unique():
        cat_df = df[df["disease_category"] == category]
        count = len(cat_df)
        weight = DISEASE_BURDEN_WEIGHTS.get(category, 1)
        banked = len(cat_df[cat_df["vault_status"] == "Banked"]) if "vault_status" in cat_df.columns else 0
        strategic_value = count * weight + banked * 2

        disease_data.append({
            "Category": category,
            "Datasets": count,
            "Weight": weight,
            "Strategic Value": strategic_value,
            "Banked": banked
        })

    disease_df = pd.DataFrame(disease_data).sort_values("Strategic Value", ascending=True)

    fig = px.bar(
        disease_df,
        y="Category",
        x="Strategic Value",
        orientation="h",
        title="Disease Burden Priority Matrix",
        color="Weight",
        color_continuous_scale="RdYlGn",
        hover_data=["Datasets", "Banked"]
    )
    fig.update_layout(height=350, margin=dict(l=20, r=20, t=40, b=20))
    st.plotly_chart(fig, use_container_width=True)


def render_team_readiness(df: pd.DataFrame):
    """Render team readiness indicators."""
    if df.empty:
        return

    st.markdown("#### Team Readiness")

    # Calculate metrics for each team
    banked = len(df[df["vault_status"] == "Banked"]) if "vault_status" in df.columns else 0
    long_read = len(df[df["seq_type"].str.contains("Long-Read", na=False)]) if "seq_type" in df.columns else 0
    gut_brain = int(df["gut_brain_relevant"].sum()) if "gut_brain_relevant" in df.columns else 0

    col1, col2, col3 = st.columns(3)

    with col1:
        wet_lab_status = "Ready" if long_read >= 5 else "Needs Data" if long_read > 0 else "Blocked"
        color = "green" if wet_lab_status == "Ready" else "orange" if wet_lab_status == "Needs Data" else "red"
        st.markdown(f"**Wet Lab**")
        st.markdown(f":{color}[{wet_lab_status}]")
        st.caption(f"{long_read} long-read datasets")

    with col2:
        comp_status = "Ready" if banked >= 10 else "Needs Data" if banked > 0 else "Blocked"
        color = "green" if comp_status == "Ready" else "orange" if comp_status == "Needs Data" else "red"
        st.markdown(f"**Computational**")
        st.markdown(f":{color}[{comp_status}]")
        st.caption(f"{banked} banked datasets")

    with col3:
        clinical_status = "Ready" if gut_brain >= 5 else "Needs Data" if gut_brain > 0 else "Blocked"
        color = "green" if clinical_status == "Ready" else "orange" if clinical_status == "Needs Data" else "red"
        st.markdown(f"**Clinical Analysis**")
        st.markdown(f":{color}[{clinical_status}]")
        st.caption(f"{gut_brain} gut-brain relevant")


def render_infrastructure_metrics(df: pd.DataFrame):
    """Render infrastructure planning metrics."""
    if df.empty:
        return

    st.markdown("#### Infrastructure Planning")

    # Calculate totals
    banked_df = df[df["vault_status"] == "Banked"] if "vault_status" in df.columns else df
    total_gb = banked_df["total_gb"].sum() if "total_gb" in banked_df.columns else 0

    col1, col2, col3 = st.columns(3)

    with col1:
        st.metric("Storage Required", f"{total_gb:.1f} Gb")
        st.caption(f"Est. S3 cost: ${total_gb * 0.023:.2f}/month")

    with col2:
        # Estimate download time at 100 Mbps
        download_hours = (total_gb * 8) / (100 * 3600) if total_gb > 0 else 0
        st.metric("Download Time", f"{download_hours:.1f} hrs")
        st.caption("At 100 Mbps connection")

    with col3:
        # Estimate processing time (rough: 1 hour per 10 Gb)
        processing_hours = total_gb / 10
        st.metric("Processing Time", f"{processing_hours:.1f} hrs")
        st.caption("Est. bioinformatics pipeline")


def render_strategic_recommendations(df: pd.DataFrame, records: list[dict]):
    """Render auto-generated strategic action items."""
    if df.empty:
        return

    st.markdown("#### Strategic Recommendations")

    recommendations = []

    banked = len(df[df["vault_status"] == "Banked"]) if "vault_status" in df.columns else 0
    pending = len(df[df["vault_status"] == "Pending Review"]) if "vault_status" in df.columns else 0
    public = len(df[df["access_type"] == "Public"]) if "access_type" in df.columns else 0
    gut_brain = int(df["gut_brain_relevant"].sum()) if "gut_brain_relevant" in df.columns else 0
    long_read = len(df[df["seq_type"].str.contains("Long-Read", na=False)]) if "seq_type" in df.columns else 0

    if banked >= 10 and public >= 5:
        recommendations.append("Initiate bulk data ingestion pipeline for high-quality public datasets")

    if pending >= 5:
        recommendations.append(f"Schedule review session for {pending} pending datasets")

    if gut_brain >= 3:
        recommendations.append(f"Prioritize {gut_brain} gut-brain relevant datasets for clinical analysis")

    if long_read >= 3:
        recommendations.append(f"Coordinate with wet lab on {long_read} long-read datasets for strain isolation")

    private = len(df[df["access_type"] == "Private (Requires Access)"]) if "access_type" in df.columns else 0
    if private >= 3:
        recommendations.append(f"Consider partnership outreach for {private} restricted-access datasets")

    if not recommendations:
        recommendations.append("Continue dataset discovery - expand search criteria for more results")

    for rec in recommendations:
        st.markdown(f"- {rec}")


def render_mission_control(df: pd.DataFrame, records: list[dict]):
    """Render the Mission Control leadership dashboard."""
    st.markdown("### Mission Control")

    # Get the selected leadership view
    current_view = st.session_state.get("leadership_view", "Banking Progress")
    st.caption(f"Holobiome Leadership Dashboard - Viewing: **{current_view}**")

    if df.empty:
        st.info("Run a search to populate the Mission Control dashboard.")
        return

    # Calculate common metrics
    total = len(df)
    banked = len(df[df["vault_status"] == "Banked"]) if "vault_status" in df.columns else 0
    pending = len(df[df["vault_status"] == "Pending Review"]) if "vault_status" in df.columns else 0
    not_suitable = len(df[df["vault_status"] == "Not Suitable"]) if "vault_status" in df.columns else 0
    gut_brain = int(df["gut_brain_relevant"].sum()) if "gut_brain_relevant" in df.columns else 0
    total_gb = df["total_gb"].sum() if "total_gb" in df.columns else 0
    public_count = len(df[df["access_type"] == "Public"]) if "access_type" in df.columns else 0
    private_count = len(df[df["access_type"] == "Private (Requires Access)"]) if "access_type" in df.columns else 0

    # ======================
    # BANKING PROGRESS VIEW
    # ======================
    if current_view == "Banking Progress":
        st.markdown("#### Banking Progress - Live Status")

        # Top metrics row
        col1, col2, col3, col4, col5 = st.columns(5)
        with col1:
            st.metric("Identified", total, help="Total datasets discovered")
        with col2:
            banked_pct = (banked / total * 100) if total > 0 else 0
            st.metric("Banked", banked, delta=f"{banked_pct:.0f}%", delta_color="normal")
        with col3:
            st.metric("Pending Review", pending)
        with col4:
            st.metric("Not Suitable", not_suitable, delta_color="inverse")
        with col5:
            st.metric("Total Volume", f"{total_gb:.1f} Gb")

        st.markdown("---")

        # Progress visualization
        col1, col2 = st.columns(2)

        with col1:
            st.markdown("##### Banking Pipeline")
            if total > 0:
                fig = go.Figure(go.Funnel(
                    y=["Discovered", "Banked", "Pending", "Not Suitable"],
                    x=[total, banked, pending, not_suitable],
                    textinfo="value+percent initial"
                ))
                fig.update_layout(height=300, margin=dict(l=20, r=20, t=20, b=20))
                st.plotly_chart(fig, use_container_width=True)

        with col2:
            st.markdown("##### Status Breakdown")
            if "vault_status" in df.columns:
                status_counts = df["vault_status"].value_counts()
                colors = {"Banked": "#2ecc71", "Pending Review": "#f39c12", "Not Suitable": "#e74c3c"}
                fig = px.pie(
                    values=status_counts.values,
                    names=status_counts.index,
                    color=status_counts.index,
                    color_discrete_map=colors
                )
                fig.update_layout(height=300, margin=dict(l=20, r=20, t=20, b=20))
                st.plotly_chart(fig, use_container_width=True)

        # Action items
        st.markdown("---")
        st.markdown("##### Action Items")
        if pending > 0:
            st.warning(f"**{pending} datasets** require manual review before banking decision")
        if banked >= 5:
            st.success(f"**{banked} datasets** ready for bulk vault ingestion")
        if banked == 0:
            st.error("No datasets meet banking criteria - consider expanding search")

    # ===========================
    # DISEASE PRIORITIZATION VIEW
    # ===========================
    elif current_view == "Disease Prioritization":
        st.markdown("#### Disease Prioritization - Strategic Focus Areas")

        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Gut-Brain Relevant", gut_brain, help="Priority research area")
        with col2:
            if "disease_category" in df.columns:
                top_cat = df["disease_category"].value_counts().head(1)
                top_name = top_cat.index[0] if not top_cat.empty else "N/A"
                top_count = top_cat.values[0] if not top_cat.empty else 0
                st.metric("Top Category", top_name, delta=f"{top_count} datasets")
        with col3:
            categories = df["disease_category"].nunique() if "disease_category" in df.columns else 0
            st.metric("Disease Categories", categories)

        st.markdown("---")

        # Disease burden chart (full width)
        render_disease_burden_chart(df)

        # Category breakdown
        if "disease_category" in df.columns:
            st.markdown("---")
            st.markdown("##### Category Details")
            cat_summary = df.groupby("disease_category").agg({
                "run_accession": "count",
                "total_gb": "sum",
                "quality_score": "mean"
            }).rename(columns={
                "run_accession": "Datasets",
                "total_gb": "Data (Gb)",
                "quality_score": "Avg Score"
            }).round(1)
            cat_summary = cat_summary.sort_values("Datasets", ascending=False)
            st.dataframe(cat_summary, use_container_width=True)

    # ===================
    # DATA QUALITY VIEW
    # ===================
    elif current_view == "Data Quality":
        st.markdown("#### Data Quality - Scoring Overview")

        col1, col2, col3, col4 = st.columns(4)
        with col1:
            avg_score = df["quality_score"].mean() if "quality_score" in df.columns else 0
            st.metric("Avg Quality Score", f"{avg_score:.0f}/100")
        with col2:
            high_quality = len(df[df["quality_score"] >= 70]) if "quality_score" in df.columns else 0
            st.metric("High Quality (70+)", high_quality)
        with col3:
            with_pubmed = len(df[df["pubmed_id"].notna() & (df["pubmed_id"] != "")]) if "pubmed_id" in df.columns else 0
            st.metric("With Publications", with_pubmed)
        with col4:
            st.metric("Total Volume", f"{total_gb:.1f} Gb")

        st.markdown("---")

        col1, col2 = st.columns(2)

        with col1:
            st.markdown("##### Quality Distribution")
            if "quality_score" in df.columns:
                fig = px.histogram(
                    df, x="quality_score", nbins=20,
                    color_discrete_sequence=["#3498db"]
                )
                fig.add_vline(x=70, line_dash="dash", line_color="green", annotation_text="Banked threshold")
                fig.add_vline(x=55, line_dash="dash", line_color="orange", annotation_text="Review threshold")
                fig.update_layout(height=300, margin=dict(l=20, r=20, t=40, b=20))
                st.plotly_chart(fig, use_container_width=True)

        with col2:
            st.markdown("##### Quality by Platform")
            if "seq_type" in df.columns and "quality_score" in df.columns:
                fig = px.box(df, x="seq_type", y="quality_score", color="seq_type")
                fig.update_layout(height=300, margin=dict(l=20, r=20, t=20, b=20), showlegend=False)
                st.plotly_chart(fig, use_container_width=True)

        # Quality factors
        st.markdown("---")
        st.markdown("##### Quality Scoring Factors")
        st.markdown("""
        | Factor | Max Points | Description |
        |--------|------------|-------------|
        | Sequencing Depth | 25 | Based on read count or throughput |
        | Read Length | 20-25 | Longer reads score higher |
        | Metadata Quality | 25 | Harmonized sample attributes |
        | Clinical Fields | 15 | Patient/clinical metadata |
        | Sample Relevance | 10 | Fecal/gut samples bonus |
        | Publication | 5 | Linked PubMed ID |
        """)

    # ====================
    # ACCESS STATUS VIEW
    # ====================
    elif current_view == "Access Status":
        st.markdown("#### Access Status - Public vs Private Data")

        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Public Datasets", public_count, help="Freely downloadable")
        with col2:
            st.metric("Private Datasets", private_count, help="Requires access request or partnership")
        with col3:
            public_pct = (public_count / total * 100) if total > 0 else 0
            st.metric("Public %", f"{public_pct:.0f}%")
        with col4:
            public_banked = len(df[(df["access_type"] == "Public") & (df["vault_status"] == "Banked")]) if "access_type" in df.columns and "vault_status" in df.columns else 0
            st.metric("Public & Bankable", public_banked)

        st.markdown("---")

        col1, col2 = st.columns(2)

        with col1:
            st.markdown("##### Access Distribution")
            if "access_type" in df.columns:
                access_counts = df["access_type"].value_counts()
                colors = {"Public": "#27ae60", "Private (Requires Access)": "#e74c3c"}
                fig = px.pie(
                    values=access_counts.values,
                    names=access_counts.index,
                    color=access_counts.index,
                    color_discrete_map=colors
                )
                fig.update_layout(height=300, margin=dict(l=20, r=20, t=20, b=20))
                st.plotly_chart(fig, use_container_width=True)

        with col2:
            st.markdown("##### Access by Status")
            if "access_type" in df.columns and "vault_status" in df.columns:
                cross_tab = pd.crosstab(df["access_type"], df["vault_status"])
                fig = px.bar(cross_tab, barmode="group")
                fig.update_layout(height=300, margin=dict(l=20, r=20, t=20, b=20))
                st.plotly_chart(fig, use_container_width=True)

        # Partnership opportunities
        st.markdown("---")
        st.markdown("##### Partnership Opportunities")
        if private_count > 0:
            private_banked_quality = df[(df["access_type"] == "Private (Requires Access)") & (df["quality_score"] >= 55)] if "access_type" in df.columns and "quality_score" in df.columns else pd.DataFrame()
            if len(private_banked_quality) > 0:
                st.info(f"**{len(private_banked_quality)} private datasets** have good quality scores - consider outreach for data sharing partnerships")
        else:
            st.success("All discovered datasets are publicly accessible")

    # ====================
    # INFRASTRUCTURE VIEW
    # ====================
    elif current_view == "Infrastructure":
        st.markdown("#### Infrastructure Planning - Storage & Compute")

        # Calculate estimates
        banked_gb = df[df["vault_status"] == "Banked"]["total_gb"].sum() if "vault_status" in df.columns and "total_gb" in df.columns else 0
        pending_gb = df[df["vault_status"] == "Pending Review"]["total_gb"].sum() if "vault_status" in df.columns and "total_gb" in df.columns else 0

        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Total Data Volume", f"{total_gb:.1f} Gb")
        with col2:
            st.metric("Banked Data", f"{banked_gb:.1f} Gb")
        with col3:
            st.metric("Pending Data", f"{pending_gb:.1f} Gb")
        with col4:
            # Estimate download time at 100 Mbps
            download_hours = (total_gb * 8) / 100 / 3600
            st.metric("Est. Download Time", f"{download_hours:.1f} hrs", help="At 100 Mbps")

        st.markdown("---")

        render_infrastructure_metrics(df)

        st.markdown("---")
        render_team_readiness(df)

    # Always show strategic recommendations at bottom
    st.markdown("---")
    render_strategic_recommendations(df, records)


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

    # Additional visualizations - Temporal and Platform
    st.markdown("---")

    col1, col2 = st.columns(2)

    with col1:
        render_temporal_trend(df, records)

    with col2:
        render_platform_breakdown(df)


def render_dataset_browser(df: pd.DataFrame, records: list[dict]):
    """Render the Dataset Browser tab with filterable table."""
    if df.empty:
        st.info("Run a search to see datasets.")
        return

    # Initialize session state for comparison
    if "selected_for_compare" not in st.session_state:
        st.session_state.selected_for_compare = []

    # Filters row 1
    col1, col2, col3, col4 = st.columns(4)

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

    with col4:
        favorites_only = st.checkbox("Show favorites only", value=False)

    # Date filtering row
    col1, col2 = st.columns(2)
    with col1:
        date_start = st.date_input("Collection date from", value=None, key="date_start")
    with col2:
        date_end = st.date_input("Collection date to", value=None, key="date_end")

    # Apply filters
    filtered_df = df.copy()
    if "vault_status" in df.columns and status_filter:
        filtered_df = filtered_df[filtered_df["vault_status"].isin(status_filter)]
    if "access_type" in df.columns and access_filter:
        filtered_df = filtered_df[filtered_df["access_type"].isin(access_filter)]
    if fecal_only and "is_fecal_sample" in df.columns:
        filtered_df = filtered_df[filtered_df["is_fecal_sample"] == True]

    # Apply favorites filter
    if favorites_only and "favorites" in st.session_state and st.session_state.favorites:
        filtered_df = filtered_df[filtered_df["run_accession"].isin(st.session_state.favorites)]

    # Apply date filter
    if date_start or date_end:
        date_filtered_accessions = []
        for record in records:
            harmonized = record.get("harmonized_fields", {})
            date_str = harmonized.get("collection_date")
            if date_str:
                parsed_date = parse_flexible_date(date_str)
                if parsed_date:
                    include = True
                    if date_start and parsed_date.date() < date_start:
                        include = False
                    if date_end and parsed_date.date() > date_end:
                        include = False
                    if include:
                        date_filtered_accessions.append(record.get("run_accession"))
        if date_filtered_accessions:
            filtered_df = filtered_df[filtered_df["run_accession"].isin(date_filtered_accessions)]

    st.caption(f"Showing {len(filtered_df)} of {len(df)} datasets")

    # Dataset comparison selector
    available_accessions = filtered_df["run_accession"].dropna().tolist() if "run_accession" in filtered_df.columns else []
    if available_accessions:
        selected_compare = st.multiselect(
            "Select datasets to compare (max 3)",
            options=available_accessions,
            default=st.session_state.selected_for_compare[:3] if st.session_state.selected_for_compare else [],
            max_selections=3,
            key="compare_selector"
        )
        st.session_state.selected_for_compare = selected_compare

    # Display table with star column for favorites
    display_df = filtered_df.copy()

    # Add favorites column
    if "favorites" not in st.session_state:
        st.session_state.favorites = set()

    if "run_accession" in display_df.columns:
        display_df["favorite"] = display_df["run_accession"].apply(
            lambda x: x in st.session_state.favorites
        )
        # Reorder columns to put favorite first
        cols = ["favorite"] + [c for c in display_df.columns if c != "favorite"]
        display_df = display_df[cols]

    edited_df = st.data_editor(
        display_df,
        use_container_width=True,
        height=400,
        column_config={
            "favorite": st.column_config.CheckboxColumn("Star", default=False, width="small"),
            "run_accession": st.column_config.TextColumn("Run ID"),
            "run_url": st.column_config.LinkColumn("Link", display_text="View"),
            "vault_status": st.column_config.TextColumn("Status"),
            "access_type": st.column_config.TextColumn("Access"),
            "disease_category": st.column_config.TextColumn("Disease"),
            "total_gb": st.column_config.NumberColumn("Gb", format="%.2f"),
            "total_score": st.column_config.NumberColumn("Score"),
        },
        disabled=[c for c in display_df.columns if c != "favorite"],
        key="dataset_browser_editor"
    )

    # Update favorites from edited dataframe
    if "favorite" in edited_df.columns and "run_accession" in edited_df.columns:
        new_favorites = set(edited_df[edited_df["favorite"] == True]["run_accession"].tolist())
        st.session_state.favorites = new_favorites

    # Export buttons
    col1, col2, col3 = st.columns(3)
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
    with col3:
        # Watchlist export
        if st.session_state.favorites:
            watchlist_df = df[df["run_accession"].isin(st.session_state.favorites)]
            st.download_button(
                "Download Watchlist (CSV)",
                watchlist_df.to_csv(index=False),
                "watchlist.csv",
                "text/csv"
            )

    # Dataset comparison view
    if st.session_state.selected_for_compare:
        st.markdown("---")
        st.markdown("### Dataset Comparison")
        render_comparison_view(st.session_state.selected_for_compare, records)


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

    col1, col2, col3, col4 = st.columns(4)

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

    with col4:
        # Watchlist export
        if "favorites" in st.session_state and st.session_state.favorites:
            watchlist_df = df[df["run_accession"].isin(st.session_state.favorites)]
            if not watchlist_df.empty:
                st.download_button(
                    "Watchlist (CSV)",
                    watchlist_df.to_csv(index=False),
                    "watchlist.csv",
                    "text/csv",
                    help="Your starred datasets"
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

    # Initialize session state
    if "results_df" not in st.session_state:
        st.session_state.results_df = pd.DataFrame()
    if "records" not in st.session_state:
        st.session_state.records = []
    if "favorites" not in st.session_state:
        st.session_state.favorites = set()
    if "selected_for_compare" not in st.session_state:
        st.session_state.selected_for_compare = []
    if "leadership_view" not in st.session_state:
        st.session_state.leadership_view = "Banking Progress"
    if "has_searched" not in st.session_state:
        st.session_state.has_searched = False

    # Read URL params for search sharing
    query_params = st.query_params
    url_query = query_params.get("q", "")
    url_role = query_params.get("role", "")
    url_preset = query_params.get("preset", "")
    url_max = query_params.get("max", "")
    url_scoring = query_params.get("scoring", "")

    # If URL has query params, skip landing page
    if url_query:
        st.session_state.has_searched = True

    # =====================
    # LANDING PAGE
    # =====================
    if not st.session_state.has_searched:
        # Centered landing page
        st.markdown("<br>", unsafe_allow_html=True)

        # Title centered
        st.markdown(
            "<h1 style='text-align: center;'>Microbiome Dataset Discovery</h1>",
            unsafe_allow_html=True
        )
        st.markdown(
            "<p style='text-align: center; color: gray;'>Search and explore microbiome sequencing datasets from NCBI SRA</p>",
            unsafe_allow_html=True
        )

        st.markdown("<br><br>", unsafe_allow_html=True)

        # Centered selection form
        col1, col2, col3 = st.columns([1, 2, 1])

        with col2:
            # Role selection
            role_options = ["Researcher", "Leadership"]
            role = st.selectbox(
                "I am a",
                role_options,
                index=0,
                key="landing_role"
            )

            st.markdown("<br>", unsafe_allow_html=True)

            if role == "Leadership":
                # Leadership sees dashboard focus options
                view_options = LEADERSHIP_VIEWS
                leadership_view = st.selectbox(
                    "and I want to see",
                    view_options,
                    index=0,
                    key="landing_view"
                )
                st.session_state.leadership_view = leadership_view
                preset = leadership_view
                query = "fecal[All Fields] AND microbiome[All Fields] AND human[Organism]"
            else:
                # Researcher sees search options
                search_options = [
                    "All Available Data",
                    "Mental Health",
                    "Pain Conditions",
                    "Digestive Health",
                    "Metabolic Health",
                    "Long-Read Sequencing",
                    "Shotgun Metagenomics",
                    "16S Amplicon",
                    "Large Cohorts",
                    "Clinical Trials",
                    "Custom Search",
                ]

                search_option = st.selectbox(
                    "and I'm looking for",
                    search_options,
                    index=0,
                    key="landing_search"
                )
                preset = search_option

                if search_option == "Custom Search":
                    query = st.text_area(
                        "Custom query",
                        value="fecal[All Fields] AND microbiome[All Fields]",
                        height=80,
                        help="NCBI Entrez syntax",
                        key="landing_custom"
                    )
                else:
                    query = RESEARCHER_SEARCHES.get(search_option, RESEARCHER_SEARCHES["All Available Data"])

            st.markdown("<br>", unsafe_allow_html=True)

            # Max results
            max_results = st.slider("Max results", 10, 100, 50, 10, key="landing_max")

            st.markdown("<br>", unsafe_allow_html=True)

            # Search button
            if st.button("Search Datasets", type="primary", use_container_width=True, key="landing_search_btn"):
                # Update URL params
                st.query_params["q"] = query
                st.query_params["role"] = role
                st.query_params["preset"] = preset
                st.query_params["max"] = str(max_results)
                st.query_params["scoring"] = "auto"

                # Show loading animation
                import random
                loading_placeholder = st.empty()
                loading_placeholder.markdown(
                    render_loading_animation(random.choice(LOADING_MESSAGES)),
                    unsafe_allow_html=True
                )

                id_list = search_sra(query, max_results)

                if id_list:
                    loading_placeholder.markdown(
                        render_loading_animation("Processing metadata..."),
                        unsafe_allow_html=True
                    )
                    records = fetch_sra_metadata(id_list)
                    scored_records = [calculate_quality_score(r, "auto") for r in records]
                    st.session_state.records = scored_records
                    st.session_state.results_df = create_results_dataframe(scored_records)
                    loading_placeholder.empty()
                    st.session_state.has_searched = True
                    st.rerun()
                else:
                    loading_placeholder.empty()
                    st.warning("No datasets found. Try adjusting your search.")

        # Footer on landing page
        st.markdown("<br><br><br>", unsafe_allow_html=True)
        st.markdown("---")
        col1, col2, col3 = st.columns(3)
        with col2:
            st.caption("Data from [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra)")

        return  # Stop here if on landing page

    # =====================
    # MAIN DATA INTERFACE
    # =====================

    # Sidebar for users who have searched
    with st.sidebar:
        # Back to start button
        if st.button("New Search", use_container_width=True):
            st.session_state.has_searched = False
            st.session_state.results_df = pd.DataFrame()
            st.session_state.records = []
            st.query_params.clear()
            st.rerun()

        st.markdown("---")

        # Mode selection - Discovery vs Analytics
        mode_options = ["Dataset Discovery", "Analytics Dashboard"]
        default_mode = "Analytics Dashboard" if url_role == "Leadership" else "Dataset Discovery"
        default_mode_idx = mode_options.index(default_mode) if default_mode in mode_options else 0

        mode = st.radio(
            "View",
            mode_options,
            index=default_mode_idx,
            horizontal=True,
            help="Switch between data discovery and leadership analytics"
        )

        st.markdown("---")

        if mode == "Analytics Dashboard":
            # Leadership analytics mode
            st.markdown("**Dashboard Focus**")
            view_options = LEADERSHIP_VIEWS
            try:
                default_view_idx = view_options.index(url_preset) if url_preset in view_options else 0
            except (ValueError, KeyError):
                default_view_idx = 0
            leadership_view = st.selectbox(
                "Show me",
                view_options,
                index=default_view_idx,
                help="Select which metrics to focus on"
            )
            st.session_state.leadership_view = leadership_view
            preset = leadership_view
            role = "Leadership"

            # Data source for analytics
            st.markdown("**Data Source**")
            query = st.text_area(
                "Search query",
                value=url_query if url_query else "fecal[All Fields] AND microbiome[All Fields] AND human[Organism]",
                height=60,
                help="Load data to analyze"
            )
        else:
            # Dataset Discovery mode (default)
            st.markdown("**Search For**")
            role = "Researcher"

            # Organized search options
            selectable_options = [
                "All Available Data",
                "Mental Health",
                "Pain Conditions",
                "Digestive Health",
                "Metabolic Health",
                "Long-Read Sequencing",
                "Shotgun Metagenomics",
                "16S Amplicon",
                "Large Cohorts",
                "Clinical Trials",
                "Custom Search",
            ]

            try:
                default_search_idx = selectable_options.index(url_preset) if url_preset in selectable_options else 0
            except (ValueError, KeyError):
                default_search_idx = 0

            search_option = st.selectbox(
                "I'm looking for",
                selectable_options,
                index=default_search_idx,
                help="Select what type of data to search for"
            )
            preset = search_option

            if search_option == "Custom Search":
                default_query = url_query if url_query else "fecal[All Fields] AND microbiome[All Fields]"
                query = st.text_area(
                    "Custom query",
                    value=default_query,
                    height=80,
                    help="NCBI Entrez syntax"
                )
            else:
                query = RESEARCHER_SEARCHES.get(search_option, RESEARCHER_SEARCHES["All Available Data"])

        default_max = int(url_max) if url_max and url_max.isdigit() else 50
        max_results = st.slider("Max results", 10, 100, default_max, 10)

        st.markdown("---")

        scoring_options = ["auto", "nanopore", "illumina"]
        try:
            default_scoring_idx = scoring_options.index(url_scoring) if url_scoring in scoring_options else 0
        except (ValueError, KeyError):
            default_scoring_idx = 0
        scoring_mode = st.radio(
            "Scoring",
            scoring_options,
            index=default_scoring_idx,
            format_func=lambda x: {"auto": "Auto", "nanopore": "Long-Read", "illumina": "Short-Read"}[x],
            horizontal=True
        )

        search_button = st.button("Update Search", type="primary", use_container_width=True)

        st.markdown("---")

        # Watchlist indicator
        if st.session_state.favorites:
            st.markdown(f"**Watchlist:** {len(st.session_state.favorites)} datasets")

        # About section
        with st.expander("About & Data Sources"):
            st.markdown("""
**Microbiome Dataset Discovery Dashboard**

A tool for discovering and evaluating microbiome sequencing datasets from public repositories.

**Data Sources:**
- [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) - Sequence Read Archive
- [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/) - Project metadata
- [PubMed](https://pubmed.ncbi.nlm.nih.gov/) - Linked publications

**How It Works:**
1. Searches NCBI using E-utilities API
2. Extracts and harmonizes metadata
3. Scores datasets for quality
4. Categorizes by disease area

**Accession Types:**
- **SRR**: Sequencing run ID
- **PRJNA**: BioProject ID
- **PMID**: PubMed article ID
            """)

    # Page title based on mode
    if mode == "Analytics Dashboard":
        st.title("Analytics Dashboard")
        st.caption("Leadership view - Genome banking operations and strategic metrics")
    else:
        st.title("Microbiome Dataset Discovery")
        st.caption("Search and explore microbiome sequencing datasets from NCBI SRA")

    # Execute search if button pressed
    if search_button:
        # Update URL params
        st.query_params["q"] = query
        st.query_params["role"] = role
        st.query_params["preset"] = preset
        st.query_params["max"] = str(max_results)
        st.query_params["scoring"] = scoring_mode

        # Show loading animation
        import random
        loading_placeholder = st.empty()
        loading_placeholder.markdown(
            render_loading_animation(random.choice(LOADING_MESSAGES)),
            unsafe_allow_html=True
        )

        id_list = search_sra(query, max_results)

        if id_list:
            loading_placeholder.markdown(
                render_loading_animation("Processing metadata..."),
                unsafe_allow_html=True
            )
            records = fetch_sra_metadata(id_list)
            scored_records = [calculate_quality_score(r, scoring_mode) for r in records]
            st.session_state.records = scored_records
            st.session_state.results_df = create_results_dataframe(scored_records)
            loading_placeholder.empty()
        else:
            loading_placeholder.empty()
            st.warning("No datasets found. Try adjusting your query.")

    df = st.session_state.results_df
    records = st.session_state.records

    # Different layouts for each mode
    if mode == "Analytics Dashboard":
        # Leadership Analytics Mode - Mission Control focused
        render_key_metrics(df)

        # Mission Control is the main content
        render_mission_control(df, records)

        # Additional tabs for details
        st.markdown("---")
        st.markdown("### Detailed Views")
        tab1, tab2, tab3, tab4 = st.tabs([
            "Disease Cohorts", "Data Quality", "Browse Data", "Export"
        ])
        with tab1:
            render_disease_cohorts(df)
        with tab2:
            render_data_quality(df)
        with tab3:
            render_dataset_browser(df, records)
        with tab4:
            render_export_tab(df, records)
    else:
        # Dataset Discovery Mode - Search & Browse focused
        render_key_metrics(df)

        # Tabs for discovery workflow
        tab1, tab2, tab3, tab4, tab5 = st.tabs([
            "Browse Datasets", "Overview", "Disease Categories", "Quality Details", "Export"
        ])
        with tab1:
            render_dataset_browser(df, records)
        with tab2:
            render_overview_tab(df, records)
        with tab3:
            render_disease_cohorts(df)
        with tab4:
            render_data_quality(df)
        with tab5:
            render_export_tab(df, records)

    # Footer
    st.markdown("---")
    st.caption(
        "Data from [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) | "
        "See 'About & Data Sources' in sidebar for more information"
    )


if __name__ == "__main__":
    main()
