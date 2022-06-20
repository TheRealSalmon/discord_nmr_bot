# bot.py
import os
import random
from tempfile import TemporaryDirectory
from dotenv import load_dotenv  # type: ignore

import discord  # type: ignore
from discord.ext import commands  # type: ignore
from rdkit import Chem  # type: ignore
from rdkit.Chem import Draw  # type: ignore

from utils import resolve_identifer


load_dotenv()
TOKEN = os.getenv('DISCORD_TOKEN')

# 2
bot = commands.Bot(command_prefix='!')

@bot.event
async def on_ready():
    print(f'{bot.user.name} has connected to Discord!')

@bot.command(name='show_mol')
async def show_mol(ctx, *identifier: str):
    if len(identifier) != 1:
        await ctx.send('only accepts one parameter')
        return
    else:
        identifier = identifier[0]
    smi = resolve_identifer(identifier)
    if smi == 'identifier could not be resolved':
        await ctx.send('invalid identifier')
        return
    mol = Chem.MolFromSmiles(smi)

    with TemporaryDirectory() as tmp:
        Chem.Draw.MolToFile(
            mol,
            f'{tmp}/mol_img.png',
            size=(300, 200),
            imageType='png'
        )
        await ctx.send(file=discord.File(f'{tmp}/mol_img.png'))

@bot.command(name='nmr_search')
async def nmr_search(ctx, *identifier: str):
    if len(identifier) != 1:
        await ctx.send('only accepts one parameter')
        return
    else:
        identifier = identifier[0]
    smi = resolve_identifer(identifier)
    if smi == 'identifier could not be resolved':
        await ctx.send('invalid identifier')
        return
    mol = Chem.MolFromSmiles(smi)
    std_inchikey = Chem.MolToInchiKey(mol)
    taut_inchikey = Chem.MolToInchiKey(mol, '/FixedH')

    await ctx.send(f'{std_inchikey}, {taut_inchikey}')

bot.run(TOKEN)
