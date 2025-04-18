import nltk
import os
from fastapi import FastAPI, Request, Form
from fastapi.responses import HTMLResponse
from Bio import Entrez, Medline
import logging
from fastapi.templating import Jinja2Templates

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

# Set up Jinja2 templates with the correct folder name
templates = Jinja2Templates(directory="Templates")

@app.get("/", response_class=HTMLResponse)
async def read_root(request: Request):
    return templates.TemplateResponse("index.html", {"request": request, "results": [], "term": "", "count": 0})

@app.get("/pubmed", response_class=HTMLResponse)
async def get_pubmed_abstracts(request: Request, term: str = "diabetes", sentences: int = 2):
    logger.info("Fetching PubMed abstracts for term: %s with %d sentences", term, sentences)
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
                            sentences_list = nltk.sent_tokenize(abstract_clean)
                            summary = " ".join(sentences_list[:sentences]) if len(sentences_list) >= sentences else "Summary unavailable"
                            abstracts.append({"title": title, "abstract": abstract, "summary": summary})
                        except Exception as e:
                            logger.error(f"Failed to process abstract for {title}: {e}")
                            raise
                    else:
                        logger.warning(f"No abstract for article: {title}")
                        abstracts.append({"title": title, "abstract": "No abstract available", "summary": "None"})

        if not abstracts:
            logger.warning("No abstracts found")
            return templates.TemplateResponse("index.html", {"request": request, "results": [{"message": "No abstracts found"}], "term": term, "count": 0})

        return templates.TemplateResponse("index.html", {"request": request, "results": abstracts, "term": term, "count": len(abstracts)})

    except Exception as e:
        logger.error(f"Error fetching PubMed data: {e}")
        return templates.TemplateResponse("index.html", {"request": request, "results": [{"error": str(e)}], "term": term, "count": 0})

@app.post("/pubmed", response_class=HTMLResponse)
async def search_pubmed(request: Request, term: str = Form(...), sentences: int = Form(2)):
    return await get_pubmed_abstracts(request, term, sentences)