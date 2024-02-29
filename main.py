
from fastapi import FastAPI, HTTPException
from dnautils import DnaUtils
from model import Dna
import json
import pprint

from database import (
    fetch_one_dna,
    fetch_all_dnas,
    create_dna,
    update_dna,
    remove_dna,
)
from fastapi.middleware.cors import CORSMiddleware

app = FastAPI()
origins = [
    "https://cute-sunflower-c9e335.netlify.app",
]
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/")
async def read_root():
    return {"Hello": "Mate"}

@app.get("/api/dna")
async def get_dna():
    response = await fetch_all_dnas()
    return response

@app.get("/api/dna/{fasta_id}", response_model=Dna)
async def get_dna_by_fasta_id(fasta_id):
    response = await fetch_one_dna(fasta_id)
    if response:
        return response
    raise HTTPException(404, f"There is no DNA with the FASTA Id {fasta_id}")

@app.post("/api/dna/", response_model=Dna)
async def post_dna(dna: Dna):
    response = await create_dna(dna.dict())
    if response:
        return response
    raise HTTPException(400, "Something went wrong")

@app.put("/api/dna/{fasta_id}/", response_model=Dna)
async def put_dna(fasta_id: str, seq: str):
    response = await update_dna(fasta_id, seq)
    if response:
        return response
    raise HTTPException(404, f"There is no DNA with the FASTA Id {fasta_id}")

@app.delete("/api/dna/{fasta_id}")
async def delete_dna(fasta_id):
    response = await remove_dna(fasta_id)
    if response:
        return "Successfully deleted DNA"
    raise HTTPException(404, f"There is no DNA with the FASTA Id {fasta_id}")


@app.post("/api/dna/utils/reverse_complement")
async def get_reverse_complement(dna: Dna):
    dna_utils = DnaUtils(dna.seq)
    reverse_complement = dna_utils.reverse_complement()
    if reverse_complement:
        return reverse_complement
    raise HTTPException(400, f"Cannot get reverse complent for sequence {dna.seq}")

@app.post("/api/dna/utils/transcript")
async def transcript(seq: str):
    dna_utils = DnaUtils(seq)
    rna = dna_utils.transcription()
    if rna:
        return rna
    raise HTTPException(400, f"Cannot transcript sequence {seq}")

@app.post("/api/dna/utils/nucleotide_freq")
async def get_nucleotide_freq(dna: Dna):
    dna_utils = DnaUtils(dna.seq)
    nuc_frq = dna_utils.nucleotide_frequemcy()
    if nuc_frq:
        return json.dumps(nuc_frq)
    raise HTTPException(400, f"Cannot process nucleotide freq for sequence {dna.seq}")

@app.post("/api/dna/utils/cg_freq")
async def get_cg_freq(dna: Dna):
    dna_utils = DnaUtils(dna.seq)
    cg_frq = dna_utils.gc_content()
    if cg_frq:
        return f'CG Content is {cg_frq} %'
    raise HTTPException(400, f"Cannot process CG freq for sequence {dna.seq}")

@app.post("/api/dna/utils/translate_to_aminoacids")
async def get_aa_seq(dna: Dna):
    dna_utils = DnaUtils(dna.seq)
    aa_seq = dna_utils.translate_to_aminoacids()
    if aa_seq:
        aas = "".join(aa_seq)
        return aas
    raise HTTPException(400, f"Cannot process Amino Acids for sequence {dna.seq}")


@app.post("/api/dna/utils/toframes")
async def get_frames(dna: Dna):
    dna_utils = DnaUtils(dna.seq)
    frames = dna_utils.gen_reading_frames()
    if frames:
        res = []
        for rf in frames:
            frame = "".join(rf)
            res.append(frame)
        return pprint.pformat(res)
    raise HTTPException(400, f"Cannot process frames for sequence {dna.seq}")

@app.post("/api/dna/utils/proteins")
async def get_proteins(seq: str):
    dna_utils = DnaUtils(seq)
    #TODO refactor
    frames = dna_utils.gen_reading_frames()
    if frames:
        res = []
        for rf in frames:
            frame = "".join(rf)
            res.append(frame)
        proteines = dna_utils.proteins_from_frame(res)
        return proteines
    raise HTTPException(400, f"Cannot process proteins for sequence {seq}")

@app.post("/api/dna/utils/hamilton_distance")
async def get_hamilton_distance(seq1: str, seq2: str):
    dna_utils = DnaUtils(seq1)
    hd = dna_utils.get_hamilton_distance_loop(seq2)
    if hd:
        return f'Hamilton distance is {hd}'
    raise HTTPException(400, f"Cannot get Hamilton distance")

@app.post("/api/dna/utils/kmer")
async def get_kmer(seq: str, kmer: str):
    dna_utils = DnaUtils(seq)
    kmer_count = dna_utils.count_kmer(kmer)
    if kmer_count:
        return f'kmer occurences for {kmer} is {kmer_count}'
    raise HTTPException(400, f"Cannot get kmer occurences for sequence {seq}")

