FROM python:3.12-slim

# --- build tools for the C++ engine ---
RUN apt-get update && apt-get install -y --no-install-recommends \
    g++ \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# --- python deps (pybind11 needed at build time, flask/gunicorn at run time) ---
COPY python_app/requirements.txt ./python_app/requirements.txt
RUN pip install --no-cache-dir pybind11 -r python_app/requirements.txt

# --- copy source ---
COPY cpp ./cpp
COPY python_app ./python_app
COPY build.sh ./build.sh

# --- compile the C++ engine into python_app/ ---
RUN bash build.sh

WORKDIR /app/python_app
EXPOSE 8080
CMD ["gunicorn", "-w", "2", "-b", "0.0.0.0:8080", "app:app"]
