FROM ubuntu:24.04

RUN apt-get update && apt-get install -y \
    g++ \
    libboost-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY index_calculus.cpp .

RUN g++ -std=c++17 -O2 index_calculus.cpp -o index_calculus

CMD ["./index_calculus"]