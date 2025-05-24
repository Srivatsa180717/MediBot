import os
import json
import faiss
import numpy as np
from sentence_transformers import SentenceTransformer
from fastapi import FastAPI, Request, Form
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates
from openai import OpenAI
from dotenv import load_dotenv
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Load API key from .env
load_dotenv()
api_key = os.getenv("MISTRAL_API_KEY")
logger.info("Loaded API key: %s", "Success" if api_key else "Not found")

# Initialize FastAPI and Jinja2
app = FastAPI()
templates = Jinja2Templates(directory="Templates")

# Load SentenceTransformer model
embedding_model = SentenceTransformer("all-MiniLM-L6-v2")

# Load FAISS index and QA pairs
logger.info("Loading FAISS index and qa_pairs.json...")
faiss_index = faiss.read_index("faiss_medquad.index")
with open("qa_pairs.json", "r", encoding="utf-8") as f:
    qa_pairs = json.load(f)
logger.info(f"Loaded {len(qa_pairs)} QAs.")

# Initialize Mistral client with error handling
try:
    client = OpenAI(
        api_key=api_key,
        base_url="https://api.mistral.ai/v1"
    )
    logger.info("Mistral client initialized successfully")
except Exception as e:
    logger.error(f"Failed to initialize Mistral client: {e}")
    client = None

@app.get("/", response_class=HTMLResponse)
async def home(request: Request):
    return templates.TemplateResponse("index.html", {
        "request": request,
        "term": "",
        "summary": "",
    })

@app.post("/qa", response_class=HTMLResponse)
async def handle_qa(request: Request, term: str = Form(...)):
    if not term.strip():
        logger.warning("Empty query received")
        return templates.TemplateResponse("index.html", {
            "request": request,
            "error": "Please enter a valid question.",
            "term": "",
            "summary": "",
        })

    logger.info(f"Processing query: {term}")
    
    try:
        # Embed the query
        query_embedding = embedding_model.encode([term]).astype("float32")
        _, indices = faiss_index.search(query_embedding, k=5)
        logger.info(f"Found {len(indices[0])} relevant documents")

        # Build context from top results
        context = ""
        for idx in indices[0]:
            if idx < len(qa_pairs):
                q = qa_pairs[idx]["question"]
                a = qa_pairs[idx]["answer"]
                context += f"Q: {q}\nA: {a}\n\n"
        
        logger.info("Context built successfully")

        # Prepare prompt
        prompt = f"You are a helpful medical assistant. Use the context below to answer the user's question.\n\n{context}\nUser: {term}\nAssistant:"

        # Call Mistral API with robust error handling
        if not client:
            logger.error("Mistral client not initialized")
            return templates.TemplateResponse("index.html", {
                "request": request,
                "term": term,
                "summary": "System error: API client not initialized properly.",
            })
            
        try:
            logger.info("Calling Mistral API...")
            response = client.chat.completions.create(
                model="mistral-small-latest",
                messages=[
                    {"role": "system", "content": "You are a helpful medical assistant."},
                    {"role": "user", "content": prompt}
                ],
                temperature=0.5,
                max_tokens=500
            )
            
            summary = response.choices[0].message.content.strip()
            logger.info("Successfully received API response")
            
            # For debugging
            logger.info(f"API response content: {summary[:100]}...")
            
            if not summary:
                logger.warning("Empty response from API")
                summary = "I couldn't generate a response. Please try a different question."
                
        except Exception as e:
            logger.error(f"API Error: {str(e)}")
            # More specific error message for debugging
            summary = f"I'm having trouble connecting to my knowledge base. Error: {str(e)}"
            
        # Return the response
        logger.info("Returning response to template")
        return templates.TemplateResponse("index.html", {
            "request": request,
            "term": term,
            "summary": summary,
        })
        
    except Exception as e:
        logger.error(f"General error: {str(e)}")
        return templates.TemplateResponse("index.html", {
            "request": request,
            "term": term,
            "summary": f"An error occurred: {str(e)}",
        })