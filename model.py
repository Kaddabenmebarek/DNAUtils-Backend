from pydantic import BaseModel

class Dna(BaseModel):
    fasta_id: str
    seq: str