from fastapi import FastAPI
from Bio import Entrez
import os

app = FastAPI()

# Set your email for PubMed API
Entrez.email = "vatsaa99@gmail.com"

@app.get("/")
async def read_root():
    return {"message": "Hello, This is MediBot! Welcome to your first web application!"}

@app.get("/pubmed")
async def get_pubmed_abstracts():
    try:
        # Search PubMed for diabetes articles
        handle = Entrez.esearch(db="pubmed", term="diabetes", retmax=10)
        record = Entrez.read(handle)
        handle.close()
        print(f"Found {len(record['IdList'])} IDs")  # Debug: Check ID count
        abstracts = []
        if record["IdList"]:
            ids = record["IdList"]
            handle = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="abstract", retmode="text")
            data = handle.read().split("\n")
            handle.close()
            print(f"Raw data split into {len(data)} chunks")  # Debug: Check split result
            # Basic preprocessing: strip whitespace
            for abstract in data:
                if abstract.strip():
                    abstracts.append(abstract.strip())
            print(f"Processed {len(abstracts)} abstracts")  # Debug: Check final count
            # Save to a local file
            if not os.path.exists("C:/temp"):
                os.makedirs("C:/temp")
            with open("C:/temp/pubmed_data.txt", "w", encoding="utf-8") as f:
                f.write("\n".join(abstracts))
                print("File written to C:/temp/pubmed_data.txt")  # Debug: Confirm write
                #You need to tranfer the txt file from temp folder to your project directory
    except Exception as e:
        print(f"Error: {e}")  # Debug: Catch any exception
    return {"count": len(abstracts), "sample": abstracts[0] if abstracts else "No abstract has been found"}