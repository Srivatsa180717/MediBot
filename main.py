import nltk
import os
from fastapi import FastAPI
from Bio import Entrez, Medline
import logging

app = FastAPI()

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# --- NLTK Setup ---
def setup_nltk():
    resources_to_download = ['punkt']
    for resource in resources_to_download:
        try:
            nltk.data.find(f'tokenizers/{resource}/english.pickle')
            print(f"NLTK '{resource}' resource already available.")
        except LookupError:
            print(f"Downloading '{resource}' resource...")
            nltk.download(resource)
        finally:
            nltk.data.path.append("C:/nltk_data")
            print("NLTK Data Path:", nltk.data.path)

    for resource in resources_to_download:
        try:
            nltk.data.load(f'tokenizers/{resource}/english.pickle')
            print(f"NLTK '{resource}' tokenizer loaded successfully.")
        except Exception as e:
            print(f"Error loading '{resource}' tokenizer: {e}")
            raise

# Run NLTK setup when the app starts
setup_nltk()

@app.get("/")
async def read_root():
    return {"message": "Hello, This is MediBot! Welcome to your first web application!"}

@app.get("/pubmed")
async def get_pubmed_abstracts(term: str = "cancer"):
    logger.info("Fetching PubMed abstracts for term: %s", term)
    Entrez.email = "vatsaa99@gmail.com"

    try:
        # Search PubMed
        search_handle = Entrez.esearch(db="pubmed", term=term, retmax=10)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        ids = search_results.get("IdList", [])
        logger.info(f"Found {len(ids)} article IDs for term: {term}")

        abstracts = []

        if ids:
            # Fetch abstracts with proper resource management
            with Entrez.efetch(db="pubmed", id=",".join(ids), rettype="medline", retmode="text") as fetch_handle:
                records = Medline.parse(fetch_handle)

                for record in records:
                    title = record.get("TI", "").strip()
                    abstract = record.get("AB", "").strip()

                    if abstract:
                        try:
                            # Preprocess abstract to remove extra newlines
                            abstract_clean = " ".join(abstract.split())

                            # Use NLTK sent_tokenize for summarization
                            sentences = nltk.sent_tokenize(abstract_clean)
                            summary = " ".join(sentences[:2]) if len(sentences) > 1 else "Summary unavailable"
                            abstracts.append({"title": title, "abstract": abstract, "summary": summary})
                        except Exception as e:
                            logger.error(f"Failed to process abstract for {title}: {e}")
                            raise  # Propagate the error for a proper response
                    else:
                        logger.warning(f"No abstract for article: {title}")
                        abstracts.append({"title": title, "abstract": "No abstract available", "summary": "None"})

        if not abstracts:
            logger.warning("No abstracts found")
            return {"count": 0, "results": [{"message": "No abstracts found"}]}

        return {
            "count": len(abstracts),
            "results": abstracts
        }

    except Exception as e:
        logger.error(f"Error fetching PubMed data: {e}")
        return {"error": str(e)}