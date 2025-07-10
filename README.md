# 🩺 MediBot – Medical Q&A Chatbot Using FAISS and Mistral AI

MediBot is a lightweight medical chatbot that answers user questions using semantic search over the **MedQuAD dataset** and AI-generated responses from **Mistral**. It uses a **FAISS index** to retrieve the most relevant Q&A context and builds prompts dynamically for high-quality answers.

---

## 📸 Demo

### 💬 Chat UI

![MediBot Chat UI](static/chat-demo.png)

### 💻 Terminal Startup

![MediBot Running](static/startup.png)

---

## 🚀 Features

- 🔍 FAISS-powered semantic search over the MedQuAD dataset
- 📚 Uses `sentence-transformers` for embedding queries and Q&As
- 💬 Simple, responsive chat interface with HTML/CSS/JS
- ⚡ FastAPI backend for performance
- 🔐 Secure API key handling via `.env`

---

## 🗂️ Project Structure

├── main.py # FastAPI backend logic
├── prepare_medquad.py # Script to generate FAISS index and JSON
├── faiss_medquad.index # Precomputed FAISS index (can be regenerated)
├── qa_pairs.json # Extracted Q&A pairs from MedQuAD
├── Templates/
│ └── index.html # Frontend chat UI
├── static/
│ ├── chat-demo.png # Screenshot of UI
│ └── startup.png # Screenshot of terminal
├── .env.example # Template for your Mistral API key
├── .gitignore # Git ignore rules
├── requirements.txt # Python dependencies
└── README.md # This file

---

## ⚙️ Setup Instructions

### 1. Clone the Repository

```bash
git clone https://github.com/yourusername/MediBot.git
cd MediBot

2. Create and Activate a Virtual Environment

python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

3. Install Dependencies

pip install -r requirements.txt

🔐 API Key Setup

Create a .env file based on .env.example and add your Mistral API key:
MISTRAL_API_KEY=your-mistral-api-key-here

▶️ Run the App

uvicorn main:app --reload
Then open: http://127.0.0.1:8000

🧪 Regenerate FAISS Index (Optional)

To rebuild the faiss_medquad.index and qa_pairs.json:
python prepare_medquad.py
Make sure the MedQuAD data is present in the correct directory format.

🧰 Technologies Used

FastAPI

FAISS

SentenceTransformers

Mistral API

HTML, CSS, JavaScript (no frontend framework)

🙋‍♂️ Author & Credits

Created by Srivatsa Adharapurapu

### Dataset Attribution

This project uses the [MedQuAD Dataset](https://github.com/abachaa/MedQuAD.git) by [abachaa], licensed under the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

Changes: Preprocessing, QA pair generation and Vector embedding using Sentence Transformers were performed for project purposes.

🛡 License

This project is licensed under the MIT License

