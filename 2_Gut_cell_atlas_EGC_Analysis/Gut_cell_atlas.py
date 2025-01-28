#!/usr/bin/env python
# coding: utf-8
---
title: "Analysis of Enteric Glia on the Gut_atlas scRNA Seq dataset"
format: 
  html:   
    default:
      body-width: 1000px
    code-fold: TRUE
date: last-modified
toc: false
df-print: paged
author: "Jay V. Patankar"
echo: fenced
warning: false
---
# ## This notebook analyzed the enteric glia signatures on the Gut cell atlas dataset of IBD patients (Oliver et al. doi.org/10.1038/s41586-024-07571-1)

# ### Importing packages and downloading the dataset

# ::: {.callout-tip}
# ## System requirements
# 
# "This is a large dataset. Recommend atleast 64GB of RAM for this analysis"
# :::

# In[1]:


# Import necessary libraries
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Enable axis ticks globally in Scanpy
sc.settings.set_figure_params(dpi=300, frameon=True, figsize=(4, 4))
sc.settings._vector_friendly = True  # Enable better vector compatibility

# Download the vasular h5ad object from the 
get_ipython().system('wget -o ./3_Extended_Pan-GI_atlas_all_lineages_18485genes_20241119.h5ad https://cellgeni.cog.sanger.ac.uk/gutcellatlas/pangi/3_Extended_Pan-GI_atlas_all_lineages_18485genes_20241119.h5ad')

# Load the h5ad file
adata = sc.read_h5ad("./3_Extended_Pan-GI_atlas_all_lineages_18485genes_20241119.h5ad")

# View the AnnData object
print(adata)


# In[2]:


# View the expression matrix dimensions (cells x genes)
display(adata.shape)

# View metadata (observations: cell-level data)
display(adata.obs.head())

# View variable (features: gene-level data)
display(adata.var.head())


# In[3]:


# Subset to only cells with metadata value 'Glia' 
 adata_glial = adata[adata.obs["level_2_annot"] == "Glia"]

# Check the unique values in the 'cell_type' column after subsetting to ensure that subsetting worked
 adata_glial.obs["level_2_annot"].unique()

# save subsetted adata object
 adata_glial.write("./glia_gut_atlas_adata.h5ad")


# In[4]:


# Subset to only cells with metadata value 'Neural' 
 adata_neuronal = adata[adata.obs["level_1_annot"] == "Neural"]

# Check the unique values in the 'cell_type' column after subsetting to ensure that subsetting worked
 adata_neuronal.obs["level_1_annot"].unique()

# save subsetted adata object
 adata_neuronal.write("./neuronal_gut_atlas_adata.h5ad")


# In[5]:


# Load the data 
adata_glial = sc.read_h5ad("./glia_gut_atlas_adata.h5ad")

print(adata_glial)


# ### Perform QC, filetering, and compute UMAP

# In[6]:


# add mitochondrial identifyier column
adata_glial.var['mt'] = adata_glial.var.index.str.startswith('MT-')


# In[7]:


adata_glial.var


# In[8]:


import pandas as pd


# In[9]:


ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"


# In[10]:


ribo_genes = pd.read_table(ribo_url, skiprows=2, header = None)
ribo_genes


# In[11]:


adata_glial.var['ribo'] = adata_glial.var_names.isin(ribo_genes[0].values)


# In[12]:


sc.pp.calculate_qc_metrics(adata_glial, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)


# In[13]:


adata_glial.obs


# In[14]:


sc.pp.filter_genes(adata_glial, min_cells=3)


# In[15]:


adata_glial.var.sort_values('n_cells_by_counts')


# In[16]:


adata_glial.obs.sort_values('n_genes_by_counts')


# In[17]:


sc.pl.violin(adata_glial, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'], 
             jitter=0.4, multi_panel=True)


# In[18]:


import numpy as np

upper_lim = np.quantile(adata_glial.obs.n_genes_by_counts.values, .98)
#upper_lim = 3000


# In[19]:


upper_lim


# In[20]:


adata_glial = adata_glial[adata_glial.obs.n_genes_by_counts < upper_lim]


# In[21]:


adata_glial = adata_glial[adata_glial.obs.pct_counts_mt < 8]


# In[22]:


adata_glial = adata_glial[adata_glial.obs.pct_counts_ribo < 15]


# In[23]:


print(adata_glial)


# In[24]:


# View the first few rows of metadata
display(adata_glial.obs.head())
display(adata_glial.obs["tissue_fraction"].value_counts())
display(adata_glial.obs["organ_groups"].value_counts())
display(adata_glial.obs["control_vs_disease"].value_counts())
display(adata_glial.obs["disease"].value_counts())


# In[25]:


# Compute and plot UMAP
sc.pp.neighbors(adata_glial)
sc.tl.umap(adata_glial)


# In[26]:


# Create the UMAP plot
sc.pl.umap(
    adata_glial,
    color="level_2_annot",
    add_outline=True,
    legend_loc="on data",
    show=False  # Do not immediately show the plot
)

# Customize the axis ticks
plt.gca().tick_params(axis="both", which="both", length=5, labelsize=10)  # Adjust size as needed
plt.xlabel("UMAP1", fontsize=12)  # Add x-axis label
plt.ylabel("UMAP2", fontsize=12)  # Add y-axis label

# Show the updated plot
plt.show()


# In[27]:


# Create the UMAP plot
sc.pl.umap(
    adata_glial,
    color="disease",
    add_outline=True,
    show=False  # Do not immediately show the plot
)

# Customize the axis ticks
plt.gca().tick_params(axis="both", which="both", length=5, labelsize=10)  # Adjust size as needed
plt.xlabel("UMAP1", fontsize=12)  # Add x-axis label
plt.ylabel("UMAP2", fontsize=12)  # Add y-axis label

# Save the plot
plt.savefig("./umap_gut_atlas_EGC_by_disease_plot.png", dpi=300, bbox_inches="tight")  # Save with high resolution


# Show the updated plot
plt.show()


# In[28]:


# Create the UMAP plot
sc.pl.umap(
    adata_glial,
    color=["SOX10","CD74", "PLP1", "NRXN1", "RET", "PHOX2B", "CCK", "NES"],
    add_outline=True,
    vmax = 8,
    show=False  # Do not immediately show the plot
)

# Customize the axis ticks
plt.gca().tick_params(axis="both", which="both", length=5, labelsize=10)  # Adjust size as needed
plt.xlabel("UMAP1", fontsize=12)  # Add x-axis label
plt.ylabel("UMAP2", fontsize=12)  # Add y-axis label

# save the fig
plt.savefig("./umap_gut_cell_atlas_EGC_genes.png", dpi=300, bbox_inches="tight")

# Show the updated plot
plt.show()


# In[ ]:


# List all metadata columns
print(adata_glial.obs.columns)


# ### Perform Glia clustering and perform clusterwise DE

# In[30]:


# Louvain clustering
sc.tl.leiden(adata_glial, resolution=0.1)
sc.pl.umap(adata_glial, color="leiden", legend_loc="on data", size=40, add_outline=True, save="Gut_cell_atlas_EGC_Clustering.png")  # Visualize clusters


# In[31]:


# Normalize the data
sc.pp.normalize_total(adata_glial, target_sum=1e4)

# Find marker genes for clusters
sc.tl.rank_genes_groups(adata_glial, groupby="leiden", method="t-test")

# Extract DE results from adata_glial
result = adata_glial.uns['rank_genes_groups']
groups = result['names'].dtype.names  # Cluster/group names


# Create an empty list to store data for all clusters
all_data = []

# Loop through each group (cluster) and extract data
for group in groups:
    group_data = pd.DataFrame({
        'cluster': group,  # Add cluster/group name
        'gene': result['names'][group],
        'pval': result['pvals'][group],
        'logfoldchange': result['logfoldchanges'][group],
        'score': result['scores'][group]
    })
    all_data.append(group_data)

# Combine all groups into a single DataFrame
de_table = pd.concat(all_data, ignore_index=True)

# Export to CSV
de_table.to_csv("de_per_cluster_EGC_gut_cell_atlas.csv", index=False)

print("Differential expression results saved to 'de_results_per_cluster.csv'")



# In[32]:


# Generate violin plots for top 20 DE genes
sc.pl.rank_genes_groups_violin(adata_glial, n_genes=20)


# In[33]:


# Heatmap
sc.pl.heatmap(adata_glial, var_names=["PLP1", "S100B", "MPZ", "SPP1", "MAL", "NRXN1", "CD74", "STAT1", "SOCS3", "CXCL9"], 
              groupby=["organ_groups", "disease"], vmax=8, dendrogram=True)

# Dotplot
sc.pl.dotplot(adata_glial, var_names=["PLP1", "S100B", "SOX10", "MAL", "MPZ", "CCK", "CD74", "DCN", "SPP1", "GFRA3", "NRXN1", "CXCL9"], 
              groupby="leiden")


# ### Add and visualize cell death pathway realted genes

# In[34]:


# load in cell death list
cell_death = pd.read_csv("cell_death_combined_Hs_scIBD.csv")
# df to list
cell_death_combined_list = cell_death["Genes"].tolist()


# In[35]:


# Extract expression data from the adata 
# Check which genes from the list are present in adata.var
genes_in_adata = [gene for gene in cell_death_combined_list if gene in adata_glial.var_names]

# Extract expression data for these genes
expression_data = adata_glial[:, genes_in_adata].X.toarray()  # Converts to dense format if sparse


# In[36]:


import numpy as np

# Calculate the mean expression across selected genes for each cell
mean_expression = np.mean(expression_data, axis=1)  # Axis 1: across genes


# In[37]:


# add mean expression as metadata 
adata_glial.obs["cell_death_signature"] = mean_expression


# In[38]:


# List all metadata columns
print(adata_glial.obs.columns)


# In[39]:


# subset to create disease specific objects
# Subset to exclude "inutero" and "preterm" from the "disease" category
adata_glia_disease = adata_glial[~adata_glial.obs["disease"].isin(["inutero", "preterm"])].copy()




# In[40]:


# Heatmap
sc.pl.heatmap(adata_glia_disease, var_names=["PLP1", "S100B", "MPZ", "SPP1", "MAL", "NRXN1", "CD74", "STAT1", "SOCS3", "CXCL9", "NES", "ART3", "MIA", "COL8A1", "NRXN3", "NRN1", "ASPA", "LGI4"], 
              groupby=["organ_groups", "sample_category"], vmax=8, dendrogram=True)

# Dotplot
sc.pl.dotplot(adata_glia_disease, var_names=["PLP1", "S100B", "SOX10", "MAL", "MPZ", "CCK", "CD74", "DCN", "SPP1", "GFRA3", "NRXN1", "CXCL9", "NES", "ART3", "MIA", "COL8A1", "NRXN3", "NRN1", "LGI4", "ASPA"], 
              groupby="disease")


# In[41]:


# Normalize the data
sc.pp.normalize_total(adata_glia_disease, target_sum=1e4)

# Find marker genes for clusters
sc.tl.rank_genes_groups(adata_glia_disease, groupby="disease", method="t-test")

# Extract DE results from adata_glial
result = adata_glia_disease.uns['rank_genes_groups']
groups = result['names'].dtype.names  # Cluster/group names


# Create an empty list to store data for all clusters
all_data = []

# Loop through each group (cluster) and extract data
for group in groups:
    group_data = pd.DataFrame({
        'cluster': group,  # Add cluster/group name
        'gene': result['names'][group],
        'pval': result['pvals'][group],
        'logfoldchange': result['logfoldchanges'][group],
        'score': result['scores'][group]
    })
    all_data.append(group_data)

# Combine all groups into a single DataFrame
de_table = pd.concat(all_data, ignore_index=True)

# Export to CSV
de_table.to_csv("de_per_disease_EGC_gut_atlas_disease.csv", index=False)

print("Differential expression results saved to 'de_results_per_cluster.csv'")



# In[42]:


# Generate violin plots for top 20 DE genes
sc.pl.rank_genes_groups_violin(adata_glia_disease, n_genes=20, save="volcano_EGC_gut_cell_atlas_DE")


# In[43]:


sc.pl.umap(adata_glial, color = ["PLP1", "CXCL9", "HLA-DRA", "cell_death_signature"], cmap="viridis", save="umap_EGCs_gut_atlas.png", add_outline=True,
         size=50, alpha=0.7, vmax=2)




# In[ ]:





# In[44]:


sc.pl.umap(adata_glial, color = ["PLP1", "CXCL9", "HLA-DRA", "cell_death_signature"], cmap="viridis", save="umap_EGCs_gut_atlas_PLP1.png", add_outline=True,
         size=50, alpha=0.7, vmax=20)


# ### Load neuronal subsets, perform filtering and QC

# In[45]:


# Load the data 
adata_neuronal = sc.read_h5ad("./neuronal_gut_atlas_adata.h5ad")

print(adata_neuronal)


# In[46]:


# View the first few rows of metadata
display(adata_neuronal.obs.head())
display(adata_neuronal.obs["tissue_fraction"].value_counts())
display(adata_neuronal.obs["organ_groups"].value_counts())
display(adata_neuronal.obs["control_vs_disease"].value_counts())
display(adata_neuronal.obs["sample_category"].value_counts())
display(adata_neuronal.obs["disease"].value_counts())
display(adata_neuronal.obs["study"].value_counts())
display(adata_neuronal.obs["donor_category"].value_counts())


# In[47]:


# add mitochondrial identifyier column
adata_neuronal.var['mt'] = adata_neuronal.var.index.str.startswith('MT-')


# In[48]:


adata_neuronal.var['ribo'] = adata_neuronal.var_names.isin(ribo_genes[0].values)


# In[49]:


sc.pp.calculate_qc_metrics(adata_neuronal, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)


# In[50]:


adata_neuronal.obs


# In[51]:


sc.pp.filter_genes(adata_neuronal, min_cells=3)


# In[52]:


adata_neuronal.var.sort_values('n_cells_by_counts')


# In[53]:


adata_neuronal.obs.sort_values('n_genes_by_counts')


# In[54]:


sc.pl.violin(adata_neuronal, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'], 
             jitter=0.4, multi_panel=True)


# In[55]:


import numpy as np

upper_lim = np.quantile(adata_neuronal.obs.n_genes_by_counts.values, .98)
#upper_lim = 3000


# In[56]:


upper_lim


# In[57]:


adata_neuronal = adata_neuronal[adata_neuronal.obs.n_genes_by_counts < upper_lim]


# In[58]:


adata_neuronal = adata_neuronal[adata_neuronal.obs.pct_counts_mt < 5]


# In[59]:


adata_neuronal = adata_neuronal[adata_neuronal.obs.pct_counts_ribo < 17.5]


# In[60]:


print(adata_neuronal)


# In[61]:


sc.pl.violin(adata_neuronal, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'], 
             jitter=0.4, multi_panel=True)


# In[62]:


# Compute and plot UMAP
sc.pp.neighbors(adata_neuronal)
sc.tl.umap(adata_neuronal)


# In[63]:


# Create the UMAP plot
sc.pl.umap(
    adata_neuronal,
    color="level_2_annot",
    add_outline=True,
    size = 30,
    show=False  # Do not immediately show the plot
)

# Customize the axis ticks
plt.gca().tick_params(axis="both", which="both", length=5, labelsize=10)  # Adjust size as needed
plt.xlabel("UMAP1", fontsize=12)  # Add x-axis label
plt.ylabel("UMAP2", fontsize=12)  # Add y-axis label

# save figure
plt.savefig("./umap_neuronal_lineages_cell_atlas.png", dpi = 300, bbox_inches="tight")

# Show the updated plot
plt.show()


# In[64]:


# Create the UMAP plot
sc.pl.umap(
    adata_neuronal,
    color="level_3_annot",
    add_outline=True,
    size = 30,
    show=False  # Do not immediately show the plot
)

# Customize the axis ticks
plt.gca().tick_params(axis="both", which="both", length=5, labelsize=10)  # Adjust size as needed
plt.xlabel("UMAP1", fontsize=12)  # Add x-axis label
plt.ylabel("UMAP2", fontsize=12)  # Add y-axis label

# save figure
plt.savefig("./umap_neuronal_lineages_cell_atlas_ontology.png", dpi = 300, bbox_inches="tight")

# Show the updated plot
plt.show()


# In[65]:


# Create the UMAP plot
sc.pl.umap(
    adata_neuronal,
    color="organ_groups",
    add_outline=True,
    size = 30,
    show=False  # Do not immediately show the plot
)

# Customize the axis ticks
plt.gca().tick_params(axis="both", which="both", length=5, labelsize=10)  # Adjust size as needed
plt.xlabel("UMAP1", fontsize=12)  # Add x-axis label
plt.ylabel("UMAP2", fontsize=12)  # Add y-axis label

# save figure
plt.savefig("./umap_neuronal_lineages_cell_atlas_organs.png", dpi = 300, bbox_inches="tight")

# Show the updated plot
plt.show()


# In[66]:


# Create the UMAP plot
sc.pl.umap(
    adata_neuronal,
    color="sample_category",
    add_outline=True,
    size = 30,
    show=False  # Do not immediately show the plot
)

# Customize the axis ticks
plt.gca().tick_params(axis="both", which="both", length=5, labelsize=10)  # Adjust size as needed
plt.xlabel("UMAP1", fontsize=12)  # Add x-axis label
plt.ylabel("UMAP2", fontsize=12)  # Add y-axis label

# save figure
plt.savefig("./umap_neuronal_lineages_cell_atlas_status.png", dpi = 300, bbox_inches="tight")

# Show the updated plot
plt.show()


# In[67]:


# Create the UMAP plot
sc.pl.umap(
    adata_neuronal,
    color=["NES", "PHOX2B", "SNAP25", "ELAVL4", "UCHL1", "HAND2", "STMN2", "PLP1", "SOX10", "HLA-DRA"],
    add_outline=True,
    vmax = 2,
    size = 18,
    show=False  # Do not immediately show the plot
)

# Customize the axis ticks
plt.gca().tick_params(axis="both", which="both", length=5, labelsize=10)  # Adjust size as needed
plt.xlabel("UMAP1", fontsize=12)  # Add x-axis label
plt.ylabel("UMAP2", fontsize=12)  # Add y-axis label

# save the fig
plt.savefig("./umap_gut_cell_atlas_neuronal_genes.png", dpi=300, bbox_inches="tight")

# Show the updated plot
plt.show()


# In[68]:


# Create the UMAP plot
sc.pl.umap(
    adata_neuronal,
    color="disease",
    add_outline=True,
    vmax = 2,
    size = 30,
    show=False  # Do not immediately show the plot
)

# Customize the axis ticks
plt.gca().tick_params(axis="both", which="both", length=5, labelsize=10)  # Adjust size as needed
plt.xlabel("UMAP1", fontsize=12)  # Add x-axis label
plt.ylabel("UMAP2", fontsize=12)  # Add y-axis label

# save the fig
plt.savefig("./umap_gut_cell_atlas_neuronal_disease.png", dpi=300, bbox_inches="tight")

# Show the updated plot
plt.show()


# In[69]:


# Extract expression data from the adata 
# Check which genes from the list are present in adata.var
genes_in_adata_n = [gene for gene in cell_death_combined_list if gene in adata_neuronal.var_names]

# Extract expression data for these genes
expression_data_n = adata_neuronal[:, genes_in_adata_n].X.toarray()  # Converts to dense format if sparse


# In[70]:


import numpy as np

# Calculate the mean expression across selected genes for each cell
mean_expression_n = np.mean(expression_data_n, axis=1)  # Axis 1: across genes


# In[71]:


# add mean expression as metadata 
adata_neuronal.obs["cell_death_signature"] = mean_expression_n


# In[72]:


# List all metadata columns
print(adata_neuronal.obs.columns)


# ### Subset to only include IBD relevant disease labels for the analysis

# In[73]:


# subset to create disease specific objects
# Subset to exclude "inutero", "preterm", and the cancer categories from the "disease" category
adata_neuronal_disease = adata_neuronal[~adata_neuronal.obs["disease"].isin(["inutero", "preterm", "cancer_colorectal", "cancer_gastric", "neighbouring_polyps", "neighbouring_cancer"])].copy()


# In[74]:


# Create the UMAP plot
sc.pl.umap(
    adata_neuronal_disease,
    color="sample_category",
    add_outline=True,
    size = 30,
    palette={"Inflamed": "tab:red", "Non_pathological": "tab:blue", "Neighbouring_inflamed": "tab:orange"},
    show=False  # Do not immediately show the plot
)

# Customize the axis ticks
plt.gca().tick_params(axis="both", which="both", length=5, labelsize=10)  # Adjust size as needed
plt.xlabel("UMAP1", fontsize=12)  # Add x-axis label
plt.ylabel("UMAP2", fontsize=12)  # Add y-axis label

# save figure
plt.savefig("./umap_neuronal_lineages_cell_atlas_status.png", dpi = 300, bbox_inches="tight")

# Show the updated plot
plt.show()


# In[75]:


# Create the UMAP plot
sc.pl.umap(
    adata_neuronal_disease,
    color="disease",
    add_outline=True,
    vmax = 2,
    size = 30,
    show=False  # Do not immediately show the plot
)

# Customize the axis ticks
plt.gca().tick_params(axis="both", which="both", length=5, labelsize=10)  # Adjust size as needed
plt.xlabel("UMAP1", fontsize=12)  # Add x-axis label
plt.ylabel("UMAP2", fontsize=12)  # Add y-axis label

# save the fig
plt.savefig("./umap_gut_cell_atlas_neuronal_disease.png", dpi=300, bbox_inches="tight")

# Show the updated plot
plt.show()


# In[76]:


# Create the UMAP plot
sc.pl.umap(
    adata_neuronal_disease,
    color="organ_groups",
    add_outline=True,
    size = 30,
    show=False  # Do not immediately show the plot
)

# Customize the axis ticks
plt.gca().tick_params(axis="both", which="both", length=5, labelsize=10)  # Adjust size as needed
plt.xlabel("UMAP1", fontsize=12)  # Add x-axis label
plt.ylabel("UMAP2", fontsize=12)  # Add y-axis label

# save figure
plt.savefig("./umap_neuronal_lineages_cell_atlas_organs.png", dpi = 300, bbox_inches="tight")

# Show the updated plot
plt.show()


# In[77]:


# Create the UMAP plot
sc.pl.umap(
    adata_neuronal_disease,
    color="level_2_annot",
    add_outline=True,
    size = 30,
    show=False  # Do not immediately show the plot
)

# Customize the axis ticks
plt.gca().tick_params(axis="both", which="both", length=5, labelsize=10)  # Adjust size as needed
plt.xlabel("UMAP1", fontsize=12)  # Add x-axis label
plt.ylabel("UMAP2", fontsize=12)  # Add y-axis label

# save figure
plt.savefig("./umap_neuronal_lineages_cell_atlas.png", dpi = 300, bbox_inches="tight")

# Show the updated plot
plt.show()


# In[78]:


sc.pl.umap(adata_neuronal, color = ["PLP1", "HLA-DRA", "CXCL9", "cell_death_signature"], cmap="viridis", save="umap_neuronal_gut_atlas.png", add_outline=True,
         size=50, alpha=0.7, vmax=0.5)


# In[79]:


sc.pl.umap(adata_neuronal_disease, color = ["PLP1", "HLA-DRA", "CXCL9", "cell_death_signature"], cmap="viridis", save="umap_neuronal_gut_atlas_disease.png", add_outline=True,
         size=50, alpha=0.7, vmax=0.5)


# In[80]:


sc.pl.umap(adata_neuronal, color = ["PLP1", "HLA-DRA", "CXCL9", "cell_death_signature"], cmap="viridis", save="umap_neuronal_gut_atlas_PLP1.png", add_outline=True,
         size=50, alpha=0.7, vmax=3.5)


# In[81]:


sc.pl.umap(adata_neuronal_disease, color = ["PLP1", "HLA-DRA", "CXCL9", "cell_death_signature"], cmap="viridis", save="umap_neuronal_gut_atlas_disease_PLP1.png", add_outline=True,
         size=50, alpha=0.7, vmax=3.5)


# In[85]:


sc.pl.umap(adata_neuronal_disease, color = ["CCL2", "CXCL10", "CD74", "HLA-DRA", "HLA-DRB1"], cmap="viridis", save="umap_neuronal_gut_atlas_disease_CCL.png", add_outline=True,
         size=50, alpha=0.7, vmax=3.5)


# In[85]:


display(adata_neuronal_disease.obs.head())


# In[83]:


# split violin plots by inflammation
sc.pl.violin(adata_neuronal_disease,
            keys=["PLP1", "CXCL9", "cell_death_signature"],
            groupby="disease",
            jitter=True,
            stripplot=True,
            rotation=90, save="violins_PLP1_CXCL9_celldeath_GCA.png")


# In[90]:


# export requirements file
get_ipython().system('pip freeze > requirements_gut_cell_atlas_EGC.txt')


# In[ ]:





# In[ ]:




