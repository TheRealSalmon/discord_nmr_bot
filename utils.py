import requests  # type: ignore
from tempfile import TemporaryDirectory

import discord  # type: ignore
from rdkit import Chem  # type: ignore


def resolve_identifer(identifier: str):
    mol = Chem.MolFromSmiles(identifier, sanitize=False)
    if mol is not None:
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol)
    return nci_chem_resolver(identifier)


def nci_chem_resolver(identifier: str, id_type: str = 'smiles'):
    try:
        nci_url = 'https://cactus.nci.nih.gov/chemical/structure/'
        r = requests.get(
            f'{nci_url}{identifier}/{id_type}'
        )
        if r.status_code != 200:
            return 'identifier could not be resolved'
        return r.content.decode()
    except:
        return 'identifier could not be resolved'

async def display_mol(ctx, mol):
    with TemporaryDirectory() as tmp:
        Chem.Draw.MolToFile(
            mol,
            f'{tmp}/mol_img.png',
            size=(300, 200),
            imageType='png'
        )
        await ctx.send(file=discord.File(f'{tmp}/mol_img.png'))
