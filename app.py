from flask import Flask, render_template, request, Response
from chemistry import compute_all
import json
from prediction import check_solubility
import pandas as pd

app = Flask(__name__)


# Take the input string from user
# Get the Property values: Polar Surface Area(PSA), Number of Rotations(Nrot), Hydrogen Bond Donor(HBD), Hydrogen Bond Acceptor(HBA)
# Render the property values, solubility and image for visualization into html
@app.route("/", methods=["GET"])
def get_properties():
    smiles_query = request.args.get("query", "")
    Nrot_value = 0
    HBD_value = 0
    HBA_value = 0
    mw_value = 0
    tpsa_value = 0

    if smiles_query != "":
        properties, props = compute_all(smiles_query)
        Nrot_value = properties["nrot"]
        HBD_value = properties["hbd"]
        HBA_value = properties["hba"]
        mw_value = properties["mw"]
        tpsa_value = properties["tpsa"]
    
        label = check_solubility("filtered_SOL.csv", props)


    # Render solubility value as well
    return render_template(
        "compounds.html",
        smiles_query=smiles_query,
        Nrot_value=Nrot_value,
        HBD_value=HBD_value,
        HBA_value=HBA_value,
        mw_value=mw_value,
        tpsa_value=tpsa_value,
        label=label
    )


# Get image of molecule
# Render image of molecule
@app.route("/molecule-image", methods=["GET"], endpoint="get_molecule_image")
def get_molecule_image():
    smiles_query = request.args.get("smiles_query", "")
    img_data = []
    if smiles_query != "":
        with open("molecules.json", "r") as f:
            molecule_dict = json.load(f)

        for key, value in molecule_dict.items():
            if value == smiles_query:
                molecule_img_name = key

        image_path = "images/" + molecule_img_name + ".png"
        with open(image_path, "rb") as img_file:
            img_data = img_file.read()
    return Response(img_data, mimetype="image/png")
