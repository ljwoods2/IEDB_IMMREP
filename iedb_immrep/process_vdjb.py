import zipfile
import airr as airrc
import wget
from iedb_immrep.utils import *


class VDJdbDataset:
    def __init__(
        self,
        file=f"vdj_db_immrep25/vdjdb-2025-02-21/vdjdb.txt",
        chains=["Chain 1", "Chain 2"],
        species="human",
        mhc_class="I",
    ):
        if os.path.exists(file) is False:
            self.df = download_vdjdb_release()
        else:
            self.df = pd.read_csv(file, sep="\t")
        self.species = species
        self.mhc_class = mhc_class

        self.df = filter_vdjdb(self.df, species, mhc_class)
        if self.species == "human":
            cdr1s, cdr2s = get_cdrs()
            self.cdr1s = cdr1s
            self.cdr2s = cdr2s
        self.chains = chains
        mapping = {"Chain 1": "A", "Chain 2": "B"}
        self.v_genes = [f"v.segm{mapping[i]}" for i in self.chains]
        self.j_genes = [f"j.segm{mapping[i]}" for i in self.chains]
        self.cdr3s = [f"cdr3{mapping[i]}" for i in self.chains]
        self.df = pair_by_complex(
            self.df,
            columns=self.v_genes + self.j_genes + self.cdr3s,
            chains=chains,
            mhc_class=mhc_class,
        )
        for chain in self.chains:
            self.df = self.df.rename(
                columns={f"cdr3{mapping[chain]}": f"CDR{mapping[chain]}3"}
            )
        if self.species == "human":
            self.cdrs = (
                ["CDRA1", "CDRA2", "CDRA3"] if "Chain 1" in self.chains else []
            )
            self.cdrs += (
                ["CDRB1", "CDRB2", "CDRB3"] if "Chain 2" in self.chains else []
            )

        else:
            self.cdrs = ["CDRA3" if "Chain 1" in self.chains else []] + [
                "CDRB3" if "Chain 2" in self.chains else []
            ]
        if self.mhc_class == "I":
            self.df["epitope"] = (
                self.df["antigen.epitope"] + "~" + self.df["mhc.a"]
            )
        else:
            self.df["epitope"] = (
                self.df["antigen.epitope"]
                + "~"
                + self.df["mhc.a"]
                + "/"
                + self.df["mhc.b"]
            )

    def filter_mhc(self):
        if self.species == "human":
            if self.mhc_class == "I":
                self.df = self.df[
                    (self.df["mhc.a"].fillna("").str.contains(":"))
                ]
            elif self.mhc_class == "II":
                self.df = self.df[
                    (self.df["mhc.b"].fillna("").str.contains(":"))
                ]
        else:
            if self.mhc_class == "I":
                self.df = self.df[~self.df["mhc.a"].isnull()]
            if self.mhc_class == "II":
                self.df = self.df[~self.df["mhc.b"].isnull()]

    def standardize_genes(self):
        self.df = correct_genes(
            self.df, self.v_genes + self.j_genes, self.species
        )

    def assign_cdrs(self):
        if self.species == "human":
            self.df = assign_cdrs(
                self.df, self.cdr1s, self.cdr2s, self.v_genes
            )
        self.df = self.df.fillna(False)
        self.df = self.df[self.df[self.cdrs].all(axis=1)]

    def run_stitchr(self):
        print(f"Running Thimble...")
        chain_mapping = {"Chain 1": "A", "Chain 2": "B"}
        rename_columns = {
            "v.segmA": "TRAV",
            "v.segmB": "TRBV",
            "j.segmA": "TRAJ",
            "j.segmB": "TRBJ",
        }
        self.df = self.df.rename(columns=rename_columns)
        thimble_outputs = run_thimble(
            self.df,
            [chain_mapping[x] for x in self.chains],
            species=self.species,
        )
        self.df["TCR_name"] = range(self.df.shape[0])
        if "Chain 1" in self.chains:
            self.df["TRA"] = self.df["TCR_name"].map(thimble_outputs["A"][0])
        if "Chain 2" in self.chains:
            self.df["TRB"] = self.df["TCR_name"].map(thimble_outputs["B"][0])

    def validate(self, df):

        if self.species == "human":
            cdra_keys = ["CDRA1", "CDRA2", "CDRA3"]
            cdrb_keys = ["CDRB1", "CDRB2", "CDRB3"]
        else:

            cdra_keys = ["CDRA3"]
            cdrb_keys = ["CDRB3"]

        val_columns = []
        if "Chain 1" in self.chains:
            validated_file = validate_file(
                df, cdr_keys=cdra_keys, seq_key="TRA"
            )
            if self.species == "human":
                val_columns += ["CDRA1_val", "CDRA2_val"]
            val_columns += ["CDRA3_val"]
        else:
            validated_file = df.copy()
        if "Chain 2" in self.chains:
            validated_file = validate_file(
                validated_file,
                cdr_keys=cdrb_keys,
                seq_key="TRB",
            )
            if self.species == "human":
                val_columns += ["CDRB1_val", "CDRB2_val"]
            val_columns += ["CDRB3_val"]
        return validated_file[validated_file[val_columns].all(axis=1)].drop(
            val_columns, axis=1
        )

    def flag_10X_studies(
        self,
        df,
        exclusion_path=os.path.join(data_path, "dat/immrep_exclusion.xlsx"),
    ):
        studies = studies_to_flag(exclusion_path, method="VDJdb")
        df["10X"] = df["reference.id"].isin(studies)
        return df

    def format_for_competition(self, output_path="vdjdb_positives.csv"):
        dedup = {}
        mapping = {"Chain 1": "A", "Chain 2": "B"}

        if self.species == "human":
            fields = {
                "A": {"TRAV", "TRAJ", "TRA", "CDRA1", "CDRA2", "CDRA3"},
                "B": {"TRBV", "TRBJ", "TRB", "CDRB1", "CDRB2", "CDRB3"},
            }
        else:

            fields = {
                "A": {"TRAV", "TRAJ", "TRA", "CDRA3"},
                "B": {"TRBV", "TRBJ", "TRB", "CDRB3"},
            }
        sequence_fields = list(
            set().union(*[fields[mapping[chain]] for chain in self.chains])
        )
        self.df["unique_tcr"] = self.df.apply(
            lambda x: "_".join([x[p] for p in self.cdrs]), axis=1
        )
        if self.mhc_class == "I":
            self.df["unique_id"] = (
                self.df["unique_tcr"]
                + "_"
                + self.df["antigen.epitope"]
                + "_"
                + self.df["mhc.a"]
            )
        else:
            self.df["unique_id"] = (
                self.df["unique_tcr"]
                + "_"
                + self.df["antigen.epitope"]
                + "_"
                + self.df["mhc.a"]
                + "_"
                + self.df["mhc.b"]
            )
        self.df["full_seq_unique"] = self.df.apply(
            lambda x: ":".join(
                [
                    (
                        x[f"TR{mapping[chain]}"]
                        if pd.notna(x[f"TR{mapping[chain]}"])
                        else "NA"
                    )
                    for chain in self.chains
                ]
            ),
            axis=1,
        )
        self.df = self.flag_10X_studies(self.df)
        for k, p in self.df.groupby("unique_id"):
            for i, (x, p) in enumerate(p.groupby(f"full_seq_unique")):
                example = p.iloc[0]
                receptor_ids = ",".join([str(p) for p in set(p["complex.id"])])
                reference_ids = ",".join(
                    [str(p) for p in set(p["reference.id"])]
                )
                just_10x = set(p["10X"].unique()).difference({True}) == set()
                for phla in example["epitope"].split(","):
                    unique_key = generate_hash()
                    dedup[unique_key] = {
                        "receptor_id": receptor_ids,
                        "references": reference_ids,
                        "epitope": phla,
                        "just_10X": just_10x,
                    }
                    dedup[unique_key].update(dict(example[sequence_fields]))
        dedup = pd.DataFrame(dedup).T
        dedup = self.validate(dedup)
        dedup["mhc"] = dedup["epitope"].map(lambda x: x.split("~")[1])
        dedup["epitope"] = dedup["epitope"].map(lambda x: x.split("~")[0])

        rename_columns = {
            "epitope": "Peptide",
            "mhc": "HLA",
            "TRAV": "Va",
            "TRBV": "Vb",
            "TRAJ": "Ja",
            "TRBJ": "Jb",
            "CDRA3": "CDR3a",
            "CDRB3": "CDR3b",
            "TRA": "TCRa",
            "TRB": "TCRb",
        }
        column_order = [
            "epitope",
            "mhc",
            "TRAV",
            "TRAJ",
            "TRA",
            "CDRA3",
            "TRBV",
            "TRBJ",
            "TRB",
            "CDRB3",
        ]
        self.final_format = dedup
        self.final_format["CDR3b_extended"] = self.final_format.apply(
            lambda x: (
                "C" + x.CDRB3 + "F"
                if pd.isna(x["TRB"])
                else x["TRB"][
                    x["TRB"].index(x["CDRB3"])
                    - 1 : x["TRB"].index(x["CDRB3"])
                    + len(x["CDRB3"])
                    + 1
                ]
            ),
            axis=1,
        )
        self.final_format["CDR3a_extended"] = self.final_format.apply(
            lambda x: (
                "C" + x.CDRA3 + "F"
                if pd.isna(x["TRA"])
                else x["TRA"][
                    x["TRA"].index(x["CDRA3"])
                    - 1 : x["TRA"].index(x["CDRA3"])
                    + len(x["CDRA3"])
                    + 1
                ]
            ),
            axis=1,
        )
        self.final_format[column_order].rename(columns=rename_columns).to_csv(
            output_path, index=False
        )

    def format_airr(self, output_path="vdjdb_positives_airr.tsv"):
        df = self.final_format.copy().reset_index()
        df = df.rename(columns={"index": "cell_id"})
        chain1_cols = ["CDRA1", "CDRA2", "CDRA3", "TRAV", "TRAJ", "TRA"]
        chain2_cols = ["CDRB1", "CDRB2", "CDRB3", "TRBV", "TRBJ", "TRB"]
        renaming = {
            p: p[:3].lower() + p[-1] + "_aa"
            for p in chain1_cols
            if p != "TRAV"
        }
        renaming.update(
            {
                p: p[:3].lower() + p[-1] + "_aa"
                for p in chain2_cols
                if p != "TRBV"
            }
        )
        renaming["TRAV"] = "v_call"
        renaming["TRBV"] = "v_call"
        renaming["TRAJ"] = "j_call"
        renaming["TRBJ"] = "j_call"
        renaming["TRA"] = "sequence_aa"
        renaming["TRB"] = "sequence_aa"
        airr = []
        for k, p in df.groupby("cell_id"):
            if "Chain 1" in self.chains:
                alpha_chains = p[chain1_cols].rename(columns=renaming)
                alpha_chains["locus"] = "TRA"
                alpha_chains["epitope"] = p.iloc[0]["epitope"]
                alpha_chains["mhc"] = p.iloc[0]["mhc"]
                airr.append(alpha_chains.assign(cell_id=k))
            if "Chain 2" in self.chains:
                beta_chains = p[chain2_cols].rename(columns=renaming)
                beta_chains["locus"] = "TRB"
                beta_chains["epitope"] = p.iloc[0]["epitope"]
                beta_chains["mhc"] = p.iloc[0]["mhc"]
                airr.append(beta_chains.assign(cell_id=k))
        airr = pd.concat(airr)
        airr["productive"] = True
        self.airr = airr
        airrc.dump_rearrangement(self.airr, output_path)


def download_vdjdb_release():
    if os.path.exists("vdjdb-2025-02-21.zip"):
        os.remove("vdjdb-2025-02-21.zip")
    filename = wget.download(
        "https://github.com/antigenomics/vdjdb-db/releases/download/pyvdjdb-2025-02-21/vdjdb-2025-02-21.zip"
    )
    print(f"Unzipping...")
    with zipfile.ZipFile(filename, "r") as zip_ref:
        zip_ref.extractall("vdj_db_immrep25")
    return pd.read_csv(f"vdj_db_immrep25/vdjdb-2025-02-21/vdjdb.txt", sep="\t")


def pair_by_complex(
    df: pd.DataFrame,
    columns,
    complex_key: str = "complex.id",
    chains=["Chain 1", "Chain 2"],
    mhc_class="I",
):
    mapping = {"Chain 1": "TRA", "Chain 2": "TRB"}
    required_fields = set([mapping[x] for x in chains])
    paired = {}
    if len(required_fields) == 1:
        for row, (k, p) in enumerate(df[df[complex_key] == 0].iterrows()):
            paired[row] = {}
            if "TRA" in required_fields:
                if p["gene"] != "TRA":
                    continue
                if mhc_class == "I":
                    alpha = dict(
                        p[
                            [
                                "v.segm",
                                "j.segm",
                                "antigen.epitope",
                                "mhc.a",
                                "cdr3",
                                "reference.id",
                            ]
                        ]
                    )
                else:
                    alpha = dict(
                        p[
                            [
                                "v.segm",
                                "j.segm",
                                "antigen.epitope",
                                "mhc.a",
                                "mhc.b",
                                "cdr3",
                                "reference.id",
                            ]
                        ]
                    )
                alpha = {
                    p + "A" if p in ["v.segm", "j.segm", "cdr3"] else p: alpha[
                        p
                    ]
                    for p in alpha
                }
                paired[row].update(alpha)
            if "TRB" in required_fields:
                print(p["gene"])
                if p["gene"] != "TRB":
                    continue
                if mhc_class == "I":
                    beta = dict(
                        p[
                            [
                                "v.segm",
                                "j.segm",
                                "antigen.epitope",
                                "mhc.a",
                                "cdr3",
                                "reference.id",
                            ]
                        ]
                    )
                else:
                    beta = dict(
                        p[
                            [
                                "v.segm",
                                "j.segm",
                                "antigen.epitope",
                                "mhc.a",
                                "mhc.b",
                                "cdr3",
                                "reference.id",
                            ]
                        ]
                    )
                beta = {
                    p + "B" if p in ["v.segm", "j.segm", "cdr3"] else p: beta[
                        p
                    ]
                    for p in beta
                }
                paired[row].update(beta)
    for k, p in df[df[complex_key] != 0].groupby(complex_key):
        paired[f"complex.id.{k}"] = {}
        if "TRA" in required_fields:
            if p[p["gene"] == "TRA"].shape[0] == 0:
                continue
            if p[p["gene"] == "TRA"].shape[0] > 1:
                continue
            if mhc_class == "I":
                alpha = dict(
                    p[p["gene"] == "TRA"].iloc[0][
                        [
                            "v.segm",
                            "j.segm",
                            "antigen.epitope",
                            "mhc.a",
                            "cdr3",
                            "reference.id",
                        ]
                    ]
                )
            else:

                alpha = dict(
                    p[p["gene"] == "TRA"].iloc[0][
                        [
                            "v.segm",
                            "j.segm",
                            "antigen.epitope",
                            "mhc.a",
                            "mhc.b",
                            "cdr3",
                            "reference.id",
                        ]
                    ]
                )
            alpha = {
                p + "A" if p in ["v.segm", "j.segm", "cdr3"] else p: alpha[p]
                for p in alpha
                if p in ["v.segm", "j.segm", "cdr3"]
            }
            paired[f"complex.id.{k}"].update(alpha)
        if "TRB" in required_fields:
            if p[p["gene"] == "TRB"].shape[0] == 0:
                continue
            if p[p["gene"] == "TRB"].shape[0] > 1:
                continue
            if mhc_class == "I":
                beta = dict(
                    p[p["gene"] == "TRB"].iloc[0][
                        [
                            "v.segm",
                            "j.segm",
                            "antigen.epitope",
                            "mhc.a",
                            "cdr3",
                            "reference.id",
                        ]
                    ]
                )
            else:
                beta = dict(
                    p[p["gene"] == "TRB"].iloc[0][
                        [
                            "v.segm",
                            "j.segm",
                            "antigen.epitope",
                            "mhc.a",
                            "mhc.b",
                            "cdr3",
                            "reference.id",
                        ]
                    ]
                )
            beta = {
                p + "B" if p in ["v.segm", "j.segm", "cdr3"] else p: beta[p]
                for p in beta
            }
            paired[f"complex.id.{k}"].update(beta)
    paired = pd.DataFrame(paired).T.fillna(False)
    paired = paired.reset_index().rename(columns={"index": "complex.id"})
    paired = paired[paired[columns].all(axis=1)]
    return paired


def correct_genes(
    df, genes=["v.segmA", "v.segmB", "j.segmA", "j.segmB"], species="human"
):
    standard_species = "homosapiens" if species == "human" else "musmusculus"
    for gene in genes:
        df[gene] = df[gene].map(
            lambda x: tt.tr.standardize(x, species=standard_species)
        )
    return df


def assign_cdrs(df, cdr1s, cdr2s, v_genes=["v.segmA", "v.segmB"]):
    for gene in v_genes:
        df[f"CDR{gene[-1]}1"] = df[gene].map(cdr1s)
        df[f"CDR{gene[-1]}2"] = df[gene].map(cdr2s)
    return df


def filter_vdjdb(df, species, mhc_class):
    if species == "human":
        return df[
            (df["mhc.class"] == f"MHC{mhc_class}")
            & (df["species"] == "HomoSapiens")
        ]
    else:
        return df[
            (df["mhc.class"] == f"MHC{mhc_class}")
            & (df["species"] == "MusMusculus")
        ]


create_leaders_constants("human")

vdjdb = VDJdbDataset(species="human", mhc_class="I")
vdjdb.filter_mhc()
vdjdb.standardize_genes()
vdjdb.assign_cdrs()
vdjdb.run_stitchr()
vdjdb.format_for_competition(
    output_path="dat/HUMAN_I/vdjdb_pos_human_I_2chain.csv"
)

vdjdb = VDJdbDataset(species="human", mhc_class="II")
vdjdb.filter_mhc()
vdjdb.standardize_genes()
vdjdb.assign_cdrs()
vdjdb.run_stitchr()
vdjdb.format_for_competition(
    output_path="dat/HUMAN_II/vdjdb_pos_human_II_2chain.csv"
)

create_leaders_constants("mouse")

vdjdb = VDJdbDataset(species="mouse", mhc_class="I")
vdjdb.filter_mhc()
vdjdb.standardize_genes()
vdjdb.assign_cdrs()
vdjdb.run_stitchr()
vdjdb.format_for_competition(
    output_path="dat/MOUSE_I/vdjdb_pos_mouse_I_2chain.csv"
)
vdjdb = VDJdbDataset(species="mouse", mhc_class="II")
vdjdb.filter_mhc()
vdjdb.standardize_genes()
vdjdb.assign_cdrs()
vdjdb.run_stitchr()
vdjdb.format_for_competition(
    output_path="dat/MOUSE_II/vdjdb_pos_mouse_II_2chain.csv"
)
