# bot.py
import os
import random
from tempfile import TemporaryDirectory
from dotenv import load_dotenv  # type: ignore

import discord  # type: ignore
from discord.ext import commands  # type: ignore
import pandas as pd  # type: ignore
from rdkit import Chem  # type: ignore
from rdkit.Chem import Draw  # type: ignore

from utils import resolve_identifer, display_mol


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

    await display_mol(ctx, mol)


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

    df = pd.read_csv('lookup.csv')

    taut_inchikey = Chem.MolToInchiKey(mol, '/FixedH')
    rows = df.loc[df.taut_inchikey == taut_inchikey]
    if len(rows) == 0:
        std_inchikey = Chem.MolToInchiKey(mol)
        rows = df.loc[df.std_inchikey == std_inchikey]
    if len(rows) == 0:
        await ctx.send(f'no entry found for {identifier}')
        return
    for row in rows.to_dict('records'):
        await ctx.send(row['name'])
        await display_mol(ctx, Chem.MolFromSmiles(row['canonical_smiles']))
        dir = f'spectra/{row["taut_inchikey"]}'
        files = os.listdir(dir)
        pngs = [png for png in files if png.endswith('.png')]
        fids = [fid for fid in files if fid.endswith('.zip')]
        for png, fid in zip(pngs, fids):
            await ctx.send(file=discord.File(f'{dir}/{png}'))
            await ctx.send(file=discord.File(f'{dir}/{fid}'))


bot.run(TOKEN)
