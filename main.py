from fastapi import FastAPI
app = FastAPI()
@app.get("/")
async def read_root():
    return {"Message":"Hello, This is MediBot! Welcome to your first web application "}