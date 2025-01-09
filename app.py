from flask import Flask, render_template, request
from chemistry import compute_all

app = Flask(__name__)


# Take the input string from user
# Get the Property values: Polar Surface Area(PSA), Number of Rotations(Nrot), Hydrogen Bond Donor(HBD), Hydrogen Bond Acceptor(HBA)
@app.route("/", methods=["GET"])
def get_properties():
    smiles_query = request.args.get("query", "")
    # smiles_query = "COC1=CC(C(C2=C3C=NC=C2)=CN(C3=O)C)=CC(OC)=C1CN4CCN(CC4)CCOCCOCCOC5=CC(C6=C(N=CS6)C)=CC=C5CNC([C@@H]7C[C@H](CN7C([C@@H](NC(C8(CC8)F)=O)C(C)(C)C)=O)O)=O"
    if smiles_query != "":
        properties = compute_all(smiles_query)

    return render_template(
        "display.html",
        smiles_query=smiles_query,
        Nrot_value=properties["nrot"],
        HBD_value=properties["hbd"],
        HBA_value=properties["hba"],
    )


# Calculate Solubility of molecule in water
# Render the property values, solubility and image for visualization into html


# @app.route("/", methods=["GET"])
# def root():
#     return render_template("base.html")
