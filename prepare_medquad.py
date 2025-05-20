import os
import json
import faiss
import numpy as np
from bs4 import BeautifulSoup
from sentence_transformers import SentenceTransformer

# Load model
model = SentenceTransformer("all-MiniLM-L6-v2")

qa_pairs = []

def extract_qa_pairs_from_xml(folder):
    for root, _, files in os.walk(folder):
        for file in files:
            if file.endswith(".xml"):
                path = os.path.join(root, file)
                try:
                    with open(path, "r", encoding="utf-8") as f:
                        soup = BeautifulSoup(f.read(), "xml")
                        qapairs = soup.find_all("QAPair")
                        for pair in qapairs:
                            q_tag = pair.find("Question")
                            a_tag = pair.find("Answer")
                            if q_tag and a_tag:
                                question = q_tag.get_text(strip=True)
                                answer = a_tag.get_text(strip=True)
                                if question and answer:
                                    qa_pairs.append({"question": question, "answer": answer})
                except Exception as e:
                    print(f" Error reading {path}: {e}")

print(" Extracting QAs from XML...")
extract_qa_pairs_from_xml("MedQuAD")
print(f"Total QA pairs found: {len(qa_pairs)}")

if not qa_pairs:
    print("No Q&A pairs found. Please check the XML file structure.")
    exit()

# Save QAs
with open("qa_pairs.json", "w", encoding="utf-8") as f:
    json.dump(qa_pairs, f, indent=2)

# Encode questions
questions = [item["question"] for item in qa_pairs]
print(" Encoding questions...")
embeddings = model.encode(questions, show_progress_bar=True).astype("float32")

# FAISS Index
index = faiss.IndexFlatL2(embeddings.shape[1])
index.add(embeddings)
faiss.write_index(index, "faiss_medquad.index")

print("All done! Saved faiss_medquad.index and qa_pairs.json")
