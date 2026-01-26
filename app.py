"""
Microbiome Dataset Discovery Dashboard

A Streamlit application for searching, scoring, and analyzing
microbiome datasets from NCBI SRA.
"""

import streamlit as st
import pandas as pd
import requests
import xml.etree.ElementTree as ET
from typing import Optional
import plotly.express as px
import plotly.graph_objects as go
from io import StringIO
import time


# Configuration
NCBI_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DEFAULT_SEARCH_TERM = "microbiome[All Fields] AND amplicon[All Fields]"
MAX_RESULTS = 100


def search_sra(query: str, max_results: int = 50) -> list[str]:
    """
    Search NCBI SRA for datasets matching the query.
    Returns a list of SRA accession IDs.
    """
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

        id_list = data.get("esearchresult", {}).get("idlist", [])
        return id_list
    except requests.RequestException as e:
        st.error(f"Search failed: {e}")
        return []


def fetch_sra_metadata(id_list: list[str]) -> list[dict]:
    """
    Fetch detailed metadata for a list of SRA IDs.
    """
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
    """
    Parse SRA XML response into structured records.
    """
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
    """
    Extract relevant fields from an EXPERIMENT_PACKAGE element.
    """
    record = {}

    # Experiment info
    experiment = exp_package.find(".//EXPERIMENT")
    if experiment is not None:
        record["accession"] = experiment.get("accession", "N/A")
        record["title"] = get_text(experiment, ".//TITLE")

        # Platform info
        platform = experiment.find(".//PLATFORM")
        if platform is not None:
            for child in platform:
                record["platform"] = child.tag
                record["instrument"] = get_text(child, "INSTRUMENT_MODEL")
                break

    # Study info
    study = exp_package.find(".//STUDY")
    if study is not None:
        record["study_accession"] = study.get("accession", "N/A")
        record["study_title"] = get_text(study, ".//STUDY_TITLE")
        record["study_abstract"] = get_text(study, ".//STUDY_ABSTRACT")

    # Sample info
    sample = exp_package.find(".//SAMPLE")
    if sample is not None:
        record["sample_accession"] = sample.get("accession", "N/A")
        record["sample_title"] = get_text(sample, ".//TITLE")
        record["organism"] = get_text(sample, ".//SCIENTIFIC_NAME")

        # Sample attributes (metadata completeness)
        attributes = sample.findall(".//SAMPLE_ATTRIBUTE")
        record["metadata_fields"] = len(attributes)
        record["sample_attributes"] = {
            get_text(attr, "TAG"): get_text(attr, "VALUE")
            for attr in attributes
        }

    # Run info (sequencing stats)
    runs = exp_package.findall(".//RUN")
    total_spots = 0
    total_bases = 0
    run_count = 0

    for run in runs:
        run_count += 1
        spots = run.get("total_spots", "0")
        bases = run.get("total_bases", "0")
        try:
            total_spots += int(spots)
            total_bases += int(bases)
        except ValueError:
            pass

    record["run_count"] = run_count
    record["total_spots"] = total_spots
    record["total_bases"] = total_bases
    record["avg_read_length"] = total_bases / total_spots if total_spots > 0 else 0

    # Library info
    library = exp_package.find(".//LIBRARY_DESCRIPTOR")
    if library is not None:
        record["library_strategy"] = get_text(library, "LIBRARY_STRATEGY")
        record["library_source"] = get_text(library, "LIBRARY_SOURCE")
        record["library_selection"] = get_text(library, "LIBRARY_SELECTION")

    return record


def get_text(element, path: str) -> str:
    """Safely extract text from an XML element."""
    if element is None:
        return "N/A"
    found = element.find(path)
    return found.text.strip() if found is not None and found.text else "N/A"


def calculate_quality_score(record: dict) -> dict:
    """
    Calculate a quality score for a dataset based on multiple factors.
    Returns the record with added quality metrics.
    """
    scores = {}

    # Sequencing depth score (0-30 points)
    # Good microbiome studies typically have >10,000 reads
    spots = record.get("total_spots", 0)
    if spots >= 100000:
        scores["depth_score"] = 30
    elif spots >= 50000:
        scores["depth_score"] = 25
    elif spots >= 10000:
        scores["depth_score"] = 20
    elif spots >= 5000:
        scores["depth_score"] = 15
    elif spots >= 1000:
        scores["depth_score"] = 10
    else:
        scores["depth_score"] = 5

    # Read length score (0-20 points)
    avg_length = record.get("avg_read_length", 0)
    if avg_length >= 250:
        scores["length_score"] = 20
    elif avg_length >= 150:
        scores["length_score"] = 15
    elif avg_length >= 100:
        scores["length_score"] = 10
    else:
        scores["length_score"] = 5

    # Metadata completeness score (0-30 points)
    metadata_fields = record.get("metadata_fields", 0)
    if metadata_fields >= 10:
        scores["metadata_score"] = 30
    elif metadata_fields >= 7:
        scores["metadata_score"] = 25
    elif metadata_fields >= 5:
        scores["metadata_score"] = 20
    elif metadata_fields >= 3:
        scores["metadata_score"] = 15
    else:
        scores["metadata_score"] = 10

    # Platform score (0-20 points)
    platform = record.get("platform", "").upper()
    if "ILLUMINA" in platform:
        scores["platform_score"] = 20
    elif "PACBIO" in platform or "OXFORD" in platform:
        scores["platform_score"] = 18
    elif "ION" in platform:
        scores["platform_score"] = 15
    else:
        scores["platform_score"] = 10

    # Calculate total score
    total_score = sum(scores.values())
    scores["total_score"] = total_score
    scores["quality_grade"] = get_quality_grade(total_score)

    record.update(scores)
    return record


def get_quality_grade(score: int) -> str:
    """Convert numeric score to letter grade."""
    if score >= 90:
        return "A"
    elif score >= 80:
        return "B"
    elif score >= 70:
        return "C"
    elif score >= 60:
        return "D"
    else:
        return "F"


def create_results_dataframe(records: list[dict]) -> pd.DataFrame:
    """
    Create a pandas DataFrame from the scored records.
    """
    if not records:
        return pd.DataFrame()

    display_columns = [
        "accession",
        "study_accession",
        "title",
        "organism",
        "platform",
        "library_strategy",
        "total_spots",
        "avg_read_length",
        "metadata_fields",
        "total_score",
        "quality_grade",
        "depth_score",
        "length_score",
        "metadata_score",
        "platform_score"
    ]

    df = pd.DataFrame(records)

    # Select only columns that exist
    available_columns = [col for col in display_columns if col in df.columns]
    df = df[available_columns]

    # Sort by quality score
    if "total_score" in df.columns:
        df = df.sort_values("total_score", ascending=False)

    return df


def render_quality_charts(df: pd.DataFrame):
    """
    Render quality distribution charts using Plotly.
    """
    if df.empty:
        st.warning("No data available for charts.")
        return

    col1, col2 = st.columns(2)

    with col1:
        # Quality grade distribution
        if "quality_grade" in df.columns:
            grade_counts = df["quality_grade"].value_counts().sort_index()
            fig = px.bar(
                x=grade_counts.index,
                y=grade_counts.values,
                labels={"x": "Quality Grade", "y": "Count"},
                title="Quality Grade Distribution",
                color=grade_counts.index,
                color_discrete_map={
                    "A": "#2ecc71",
                    "B": "#3498db",
                    "C": "#f39c12",
                    "D": "#e67e22",
                    "F": "#e74c3c"
                }
            )
            fig.update_layout(showlegend=False)
            st.plotly_chart(fig, use_container_width=True)

    with col2:
        # Score component breakdown
        score_cols = ["depth_score", "length_score", "metadata_score", "platform_score"]
        available_scores = [col for col in score_cols if col in df.columns]

        if available_scores:
            avg_scores = df[available_scores].mean()
            fig = px.bar(
                x=[col.replace("_score", "").title() for col in avg_scores.index],
                y=avg_scores.values,
                labels={"x": "Score Component", "y": "Average Score"},
                title="Average Score by Component",
                color=avg_scores.values,
                color_continuous_scale="Viridis"
            )
            fig.update_layout(showlegend=False)
            st.plotly_chart(fig, use_container_width=True)

    # Sequencing depth histogram
    if "total_spots" in df.columns:
        fig = px.histogram(
            df,
            x="total_spots",
            nbins=20,
            title="Sequencing Depth Distribution",
            labels={"total_spots": "Total Reads", "count": "Number of Datasets"}
        )
        fig.update_layout(showlegend=False)
        st.plotly_chart(fig, use_container_width=True)

    # Platform breakdown
    if "platform" in df.columns:
        col1, col2 = st.columns(2)

        with col1:
            platform_counts = df["platform"].value_counts()
            fig = px.pie(
                values=platform_counts.values,
                names=platform_counts.index,
                title="Sequencing Platforms"
            )
            st.plotly_chart(fig, use_container_width=True)

        with col2:
            if "library_strategy" in df.columns:
                strategy_counts = df["library_strategy"].value_counts()
                fig = px.pie(
                    values=strategy_counts.values,
                    names=strategy_counts.index,
                    title="Library Strategies"
                )
                st.plotly_chart(fig, use_container_width=True)


def convert_df_to_csv(df: pd.DataFrame) -> str:
    """Convert DataFrame to CSV string for download."""
    return df.to_csv(index=False)


def main():
    """Main application entry point."""
    st.set_page_config(
        page_title="Microbiome Dataset Discovery",
        page_icon=None,
        layout="wide"
    )

    st.title("Microbiome Dataset Discovery Dashboard")
    st.markdown(
        "Search and evaluate microbiome datasets from NCBI SRA based on "
        "sequencing quality, depth, and metadata completeness."
    )

    # Sidebar controls
    with st.sidebar:
        st.header("Search Parameters")

        search_term = st.text_input(
            "Search Query",
            value=DEFAULT_SEARCH_TERM,
            help="NCBI search query. Use field tags like [All Fields], [Organism], etc."
        )

        max_results = st.slider(
            "Maximum Results",
            min_value=10,
            max_value=MAX_RESULTS,
            value=25,
            step=5
        )

        st.markdown("---")
        st.markdown("### Example Queries")
        st.markdown("""
        - `microbiome[All Fields] AND human[Organism]`
        - `16S[All Fields] AND gut[All Fields]`
        - `metagenome[All Fields] AND soil[All Fields]`
        - `amplicon[All Fields] AND mouse[Organism]`
        """)

        search_button = st.button("Search SRA", type="primary", use_container_width=True)

    # Initialize session state
    if "results_df" not in st.session_state:
        st.session_state.results_df = pd.DataFrame()

    # Execute search
    if search_button:
        with st.spinner("Searching NCBI SRA..."):
            id_list = search_sra(search_term, max_results)

        if id_list:
            st.info(f"Found {len(id_list)} datasets. Fetching metadata...")

            with st.spinner("Fetching and scoring datasets..."):
                records = fetch_sra_metadata(id_list)

                # Score each record
                scored_records = [calculate_quality_score(r) for r in records]

                # Create DataFrame
                st.session_state.results_df = create_results_dataframe(scored_records)

            st.success(f"Successfully processed {len(st.session_state.results_df)} datasets.")
        else:
            st.warning("No results found. Try a different search query.")

    # Display results
    df = st.session_state.results_df

    if not df.empty:
        st.markdown("---")

        # Summary metrics
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Total Datasets", len(df))
        with col2:
            avg_score = df["total_score"].mean() if "total_score" in df.columns else 0
            st.metric("Avg Quality Score", f"{avg_score:.1f}")
        with col3:
            high_quality = len(df[df["quality_grade"].isin(["A", "B"])]) if "quality_grade" in df.columns else 0
            st.metric("High Quality (A/B)", high_quality)
        with col4:
            avg_depth = df["total_spots"].mean() if "total_spots" in df.columns else 0
            st.metric("Avg Read Count", f"{avg_depth:,.0f}")

        st.markdown("---")

        # Tabs for different views
        tab1, tab2 = st.tabs(["Data Table", "Quality Charts"])

        with tab1:
            st.subheader("Dataset Results")

            # Filters
            col1, col2 = st.columns(2)
            with col1:
                if "quality_grade" in df.columns:
                    grade_filter = st.multiselect(
                        "Filter by Quality Grade",
                        options=sorted(df["quality_grade"].unique()),
                        default=sorted(df["quality_grade"].unique())
                    )
            with col2:
                if "platform" in df.columns:
                    platform_filter = st.multiselect(
                        "Filter by Platform",
                        options=sorted(df["platform"].dropna().unique()),
                        default=sorted(df["platform"].dropna().unique())
                    )

            # Apply filters
            filtered_df = df.copy()
            if "quality_grade" in df.columns and grade_filter:
                filtered_df = filtered_df[filtered_df["quality_grade"].isin(grade_filter)]
            if "platform" in df.columns and platform_filter:
                filtered_df = filtered_df[filtered_df["platform"].isin(platform_filter)]

            # Display table
            st.dataframe(
                filtered_df,
                use_container_width=True,
                height=400
            )

            # Download button
            csv_data = convert_df_to_csv(filtered_df)
            st.download_button(
                label="Download CSV",
                data=csv_data,
                file_name="microbiome_datasets.csv",
                mime="text/csv"
            )

        with tab2:
            st.subheader("Quality Distribution Analysis")
            render_quality_charts(df)

    else:
        st.info(
            "Enter a search query and click 'Search SRA' to discover microbiome datasets."
        )

    # Footer
    st.markdown("---")
    st.markdown(
        "*Data sourced from NCBI Sequence Read Archive (SRA). "
        "Quality scores are calculated based on sequencing depth, read length, "
        "metadata completeness, and platform.*"
    )


if __name__ == "__main__":
    main()
