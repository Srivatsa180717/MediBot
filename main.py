from fastapi import FastAPI
from Bio import Entrez, Medline
import os
import logging

app = FastAPI()

# Logging setup
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@app.get("/")
async def read_root():
    return {"message": "Hello, This is MediBot! Welcome to your first web application!"}

@app.get("/pubmed")
async def get_pubmed_abstracts():
    Entrez.email = "vatsaa99@gmail.com"  # Required for NCBI API access

    try:
        # Search PubMed for relevant articles
        search_handle = Entrez.esearch(db="pubmed", term="diabetes", retmax=10)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        ids = search_results.get("IdList", [])
        logger.info(f"Found {len(ids)} article IDs")

        abstracts = []

        if ids:
            # Fetch full records in MEDLINE format (structured and clean)
            fetch_handle = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="medline", retmode="text")
            records = Medline.parse(fetch_handle)

            for record in records:
                title = record.get("TI", "").strip()
                abstract = record.get("AB", "").strip()
                if abstract:
                    full_text = f"{title}\n\n{abstract}"
                    abstracts.append(full_text)
            fetch_handle.close()

        logger.info(f"Processed {len(abstracts)} abstracts")

        # Save to file
        output_dir = "C:/temp"
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, "pubmed_data.txt")
        with open(output_path, "w", encoding="utf-8") as f:
            f.write("\n\n---\n\n".join(abstracts))
        logger.info(f"Saved abstracts to {output_path}")

        return {
            "count": len(abstracts),
            "sample": abstracts[0] if abstracts else "No abstract found"
        }

    except Exception as e:
        logger.error(f"Error fetching PubMed data: {e}")
        return {"error": str(e)}
