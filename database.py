import motor.motor_asyncio
from model import Dna


uri = "mongodb+srv://dnautilz:Roche!Bobois41*@cluster0.qfert5g.mongodb.net/?retryWrites=true&w=majority&appName=Cluster0"
client = motor.motor_asyncio.AsyncIOMotorClient(uri)
#client = motor.motor_asyncio.AsyncIOMotorClient('mongodb://localhost:27017/')
database = client.DNAUtils
collection = database.DNA

async def fetch_one_dna(fasta_id):
    document = await collection.find_one({"fasta_id": fasta_id})
    return document

async def fetch_all_dnas():
    dnas = []
    cursor = collection.find({})
    async for document in cursor:
        dnas.append(Dna(**document))
    return dnas

async def create_dna(dna):
    document = dna
    result = await collection.insert_one(document)
    return document


async def update_dna(fasta_id, seq):
    await collection.update_one({"fasta_id": fasta_id}, {"$set": {"seq": seq}})
    document = await collection.find_one({"fasta_id": fasta_id})
    return document

async def remove_dna(fasta_id):
    await collection.delete_one({"fasta_id": fasta_id})
    return True