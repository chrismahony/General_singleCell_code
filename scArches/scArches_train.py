import numpy as np
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import classification_report
from scarches.models.scpoli import scPoli

import warnings
warnings.filterwarnings('ignore')
%load_ext autoreload
%autoreload 2
%matplotlib inline


  !cd ~
!ln -s /rds/projects/c/croftap-celldive01/amp2/prepare_scArches scArches_amp2_mapjag

  fibs_amp2 = sc.read('/rds/projects/c/croftap-celldive01/amp2/prepare_scArches/fibs_amp2_seurat.h5ad')
  fibs_MPAJAG = sc.read('/rds/projects/c/croftap-celldive01/amp2/prepare_scArches/fibs_mapjag_seurat.h5ad')

sc.pl.umap(fibs_MPAJAG, color=['clusters'], wspace=0.5, frameon=False)
sc.pl.umap(fibs_amp2, color=['cluster_name'], wspace=0.5, frameon=False)

early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
} 

sc.pp.neighbors(fibs_MPAJAG)

cell_type_key='clusters'
condition_key='orig.ident'

scpoli_model = scPoli(
    adata=fibs_MPAJAG,
    condition_keys=condition_key,
    cell_type_keys=cell_type_key,
    embedding_dims=10,
    recon_loss='mse',
)

scpoli_model.train(
    n_epochs=50,
    pretraining_epochs=40,
    early_stopping_kwargs=early_stopping_kwargs,
    eta=5,
)    


fibs_amp2.obs['clusters']=fibs_amp2.obs['cluster_name']
fibs_amp2.obs['orig.ident']=fibs_amp2.obs['sample']
fibs_amp2

scpoli_query = scPoli.load_query_data(
    adata=fibs_amp2,
    reference_model=scpoli_model,
    labeled_indices=[],
    )

scpoli_query.train(
    n_epochs=50,
    pretraining_epochs=40,
    eta=10
)

fibs_amp2.X = fibs_amp2.X.astype('float32')
results_dict = scpoli_query.classify(fibs_amp2, scale_uncertainties=True)

for i in range(len(cell_type_key)):
    preds = results_dict[cell_type_key]["preds"]
    results_dict[cell_type_key]["uncert"]
    classification_df = pd.DataFrame(
        classification_report(
            y_true=fibs_amp2.obs[cell_type_key],
            y_pred=preds,
            output_dict=True,
        )
    ).transpose()
print(classification_df)


scpoli_query.model.eval()
data_latent_source = scpoli_query.get_latent(
    fibs_MPAJAG,
    mean=True
)

adata_latent_source = sc.AnnData(data_latent_source)
adata_latent_source.obs = fibs_MPAJAG.obs.copy()

data_latent= scpoli_query.get_latent(
    fibs_amp2,
    mean=True
)

adata_latent = sc.AnnData(data_latent)
adata_latent.obs = fibs_amp2.obs.copy()


adata_latent.obs['cell_type_pred'] = results_dict['clusters']['preds'].tolist()
adata_latent.obs['cell_type_uncert'] = results_dict['clusters']['uncert'].tolist()
adata_latent.obs['classifier_outcome'] = (
    adata_latent.obs['cell_type_pred'] == adata_latent.obs['clusters']
)

labeled_prototypes = scpoli_query.get_prototypes_info()
labeled_prototypes.obs['study'] = 'labeled prototype'
unlabeled_prototypes = scpoli_query.get_prototypes_info(prototype_set='unlabeled')
unlabeled_prototypes.obs['study'] = 'unlabeled prototype'

adata_latent_full = adata_latent_source.concatenate(
    [adata_latent, labeled_prototypes, unlabeled_prototypes],
    batch_key='query'
)

adata_latent_full.obs['cell_type_pred'][adata_latent_full.obs['query'].isin(['0'])] = np.nan
sc.pp.neighbors(adata_latent_full, n_neighbors=15)
sc.tl.umap(adata_latent_full)

adata_no_prototypes = adata_latent_full[adata_latent_full.obs['query'].isin(['0', '1'])]

sc.pl.umap(
    adata_no_prototypes,
    color='cell_type_pred',
    show=False,
    frameon=False,
)                  

meta_adata_no_prototypes=adata_no_prototypes.obs
meta_adata_no_prototypes_df = pd.DataFrame(meta_adata_no_prototypes) 
meta_adata_no_prototypes_df.to_csv('/rds/projects/c/croftap-celldive01/amp2/prepare_scArches/meta_amp2_scarhces_labels_df.csv')

import dill
dill.dump_session('/rds/projects/c/croftap-celldive01/amp2/prepare_scArches/analysis_env.db')             
