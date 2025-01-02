#!/usr/bin/env python
# coding: utf-8

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


# In[ ]:


# add mitochondrial identifyier column
adata.var['mt'] = adata.var.index.str.startswith('MT-')


# In[ ]:


# import list of ribosomal genes
ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"


# In[ ]:


# load the ribosomal genes
ribo_genes = pd.read_table(ribo_url, skiprows=2, header = None)
ribo_genes


# In[ ]:


# add metadata identifyier for ribosomal genes to the object
adata.var['ribo'] = adata.var_names.isin(ribo_genes[0].values)


# In[ ]:


# calculate quality metrics such as mt and ribo percentage rna 
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)


# In[ ]:


# view the qc metadata
adata.obs


# In[ ]:


# filter for genes expressed in minimum of n cells
sc.pp.filter_genes(adata_glial, min_cells=3)


# In[ ]:


# sort expression data by number of cells by total counts
adata.var.sort_values('n_cells_by_counts')


# In[ ]:


# sort cell metadata by number of cells by total counts
adata.obs.sort_values('n_genes_by_counts')


# In[ ]:


# plot QC metrics
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'], 
             jitter=0.4, multi_panel=True)


# In[ ]:


# set 98th percentile filtering threshold for genes by counts
import numpy as np

upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)


# In[ ]:


# filter by the upper limit 
adata= adata[adata.obs.n_genes_by_counts < upper_lim]


# In[ ]:


# filter by mitochondrial threshold
adata = adata[adata.obs.pct_counts_mt < 8]


# In[ ]:


# filter by ribosomal threshold
adata = adata[adata.obs.pct_counts_ribo < 17.5]


# In[ ]:


# view the filtering changes 
print(adata)


# In[ ]:


# plot qc metrices post filtering
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'], 
             jitter=0.4, multi_panel=True)


# In[ ]:





# In[ ]:





# In[3]:


# Subset to only cells with metadata value 'Glia' 
adata_glial = adata[adata.obs["level_2_annot"] == "Glia"]

# Check the unique values in the 'cell_type' column after subsetting to ensure that subsetting worked
adata_glial.obs["level_2_annot"].unique()

# save subsetted adata object so as to work with the smaller object in the future
adata_glial.write("./glia_gut_atlas_adata.h5ad")


# In[4]:


# Subset to only cells with metadata value 'Neural' 
adata_neuronal = adata[adata.obs["level_1_annot"] == "Neural"]

# Check the unique values in the 'cell_type' column after subsetting to ensure that subsetting worked
adata_neuronal.obs["level_1_annot"].unique()

# save subsetted adata object so as to work with the smaller object in the future
adata_neuronal.write("./neuronal_gut_atlas_adata.h5ad")


# In[ ]:


# View the first few rows of metadata
display(adata_glial.obs.head())
display(adata_glial.obs["level_1_annot"].value_counts())
display(adata_glial.obs["level_2_annot"].value_counts())



# In[5]:


# check the data 
print(adata_glial)


# In[6]:


# View the first few rows of metadata
display(adata_glial.obs.head())
display(adata_glial.obs["tissue_fraction"].value_counts())
display(adata_glial.obs["organ_groups"].value_counts())
display(adata_glial.obs["control_vs_disease"].value_counts())
display(adata_glial.obs["disease"].value_counts())


# In[7]:


# Compute and plot UMAP
sc.pp.neighbors(adata_glial)
sc.tl.umap(adata_glial)


# In[8]:


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


# In[9]:


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


# In[10]:


# Create the UMAP plot
sc.pl.umap(
    adata_glial,
    color=["SOX10","CD74", "PLP1", "NRXN1", "RET", "PHOX2B", "CCK", "NES"],
    add_outline=True,
    vmax = 8,
    show=False  # Do not immediately display plot
)

# Customize the axis ticks
plt.gca().tick_params(axis="both", which="both", length=5, labelsize=10)  
plt.xlabel("UMAP1", fontsize=12)  
plt.ylabel("UMAP2", fontsize=12)

# save the fig
plt.savefig("./umap_gut_cell_atlas_EGC_genes.png", dpi=300, bbox_inches="tight")

# Show the updated plot
plt.show()


# In[12]:


# Louvain clustering
sc.tl.leiden(adata_glial, resolution=0.1)
sc.pl.umap(adata_glial, color="leiden", legend_loc="on data", size=40, add_outline=True, save="Gut_cell_atlas_EGC_Clustering.png")  # Visualize clusters


# In[13]:


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



# In[14]:


# Generate violin plots for top 20 DE genes
sc.pl.rank_genes_groups_violin(adata_glial, n_genes=20)


# In[16]:


# load in cell death list
cell_death = pd.read_csv("cell_death_combined_Hs_scIBD.csv")
# df to list
cell_death_combined_list = cell_death["Genes"].tolist()


# In[17]:


# Extract expression data from the adata 
# Check which genes from the list are present in adata.var
genes_in_adata = [gene for gene in cell_death_combined_list if gene in adata_glial.var_names]

# Extract expression data for these genes
expression_data = adata_glial[:, genes_in_adata].X.toarray()  # Converts to dense format if sparse


# In[18]:


import numpy as np

# Calculate the mean expression across selected genes for each cell
mean_expression = np.mean(expression_data, axis=1)  # Axis 1: across genes


# In[19]:


# add mean expression as metadata 
adata_glial.obs["cell_death_signature"] = mean_expression


# In[20]:


# List all metadata columns
print(adata_glial.obs.columns)


# In[40]:


# subset to create disease specific objects
# Subset to exclude "inutero" and "preterm" from the "disease" category
adata_glia_disease = adata_glial[~adata_glial.obs["disease"].isin(["inutero", "preterm"])].copy()




# In[60]:


# Heatmap
sc.pl.heatmap(adata_glia_disease, var_names=["PLP1", "S100B", "MPZ", "SPP1", "MAL", "NRXN1", "CD74", "STAT1", "SOCS3", "CXCL9", "NES", "ART3", "MIA", "COL8A1", "NRXN3", "NRN1", "ASPA", "LGI4"], 
              groupby=["organ_groups", "sample_category"], vmax=8, dendrogram=True)

# Dotplot
sc.pl.dotplot(adata_glia_disease, var_names=["PLP1", "S100B", "SOX10", "MAL", "MPZ", "CCK", "CD74", "DCN", "SPP1", "GFRA3", "NRXN1", "CXCL9", "NES", "ART3", "MIA", "COL8A1", "NRXN3", "NRN1", "LGI4", "ASPA"], 
              groupby="disease")


# In[61]:


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



# In[65]:


# Generate violin plots for top 20 DE genes
sc.pl.rank_genes_groups_violin(adata_glia_disease, n_genes=20, save="volcano_EGC_gut_cell_atlas_DE")


# In[22]:


sc.pl.umap(adata_glial, color = ["PLP1", "CXCL9", "HLA-DRA", "cell_death_signature"], cmap="viridis", save="umap_EGCs_gut_atlas.png", add_outline=True,
         size=50, alpha=0.7, vmax=2)




# In[23]:


sc.pl.umap(adata_glial, color = ["PLP1", "CXCL9", "HLA-DRA", "cell_death_signature"], cmap="viridis", save="umap_EGCs_gut_atlas_PLP1.png", add_outline=True,
         size=50, alpha=0.7, vmax=12)


# In[24]:


# Load the data 
print(adata_neuronal)


# In[25]:


# View the first few rows of metadata
display(adata_neuronal.obs.head())
display(adata_neuronal.obs["tissue_fraction"].value_counts())
display(adata_neuronal.obs["organ_groups"].value_counts())
display(adata_neuronal.obs["control_vs_disease"].value_counts())
display(adata_neuronal.obs["sample_category"].value_counts())
display(adata_neuronal.obs["disease"].value_counts())
display(adata_neuronal.obs["study"].value_counts())
display(adata_neuronal.obs["donor_category"].value_counts())


# In[26]:


# Compute and plot UMAP
sc.pp.neighbors(adata_neuronal)
sc.tl.umap(adata_neuronal)


# In[27]:


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


# In[28]:


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


# In[29]:


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


# In[30]:


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


# In[31]:


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


# In[32]:


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


# In[33]:


# Extract expression data from the adata 
# Check which genes from the list are present in adata.var
genes_in_adata_n = [gene for gene in cell_death_combined_list if gene in adata_neuronal.var_names]

# Extract expression data for these genes
expression_data_n = adata_neuronal[:, genes_in_adata_n].X.toarray()  # Converts to dense format if sparse


# In[34]:


import numpy as np

# Calculate the mean expression across selected genes for each cell
mean_expression_n = np.mean(expression_data_n, axis=1)  # Axis 1: across genes


# In[35]:


# add mean expression as metadata 
adata_neuronal.obs["cell_death_signature"] = mean_expression_n


# In[36]:


# List all metadata columns
print(adata_neuronal.obs.columns)


# In[44]:


# subset to create disease specific objects
# Subset to exclude "inutero" and "preterm" from the "disease" category
adata_neuronal_disease = adata_neuronal[~adata_neuronal.obs["disease"].isin(["inutero", "preterm"])].copy()




# In[ ]:


# Create the UMAP plot
sc.pl.umap(
    adata_neuronal_disease,
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


# In[ ]:


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


# In[ ]:


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


# In[ ]:


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


# In[49]:


sc.pl.umap(adata_neuronal, color = ["PLP1", ""HLA-DRA", "CXCL9", "cell_death_signature"], cmap="viridis", save="umap_neuronal_gut_atlas.png", add_outline=True,
         size=50, alpha=0.7, vmax=0.5)


# In[50]:


sc.pl.umap(adata_neuronal_disease, color = ["PLP1", ""HLA-DRA", "CXCL9", "cell_death_signature"], cmap="viridis", save="umap_neuronal_gut_atlas_disease.png", add_outline=True,
         size=50, alpha=0.7, vmax=0.5)


# In[51]:


sc.pl.umap(adata_neuronal, color = ["PLP1", ""HLA-DRA", "CXCL9", "cell_death_signature"], cmap="viridis", save="umap_neuronal_gut_atlas_PLP1.png", add_outline=True,
         size=50, alpha=0.7, vmax=3.5)


# In[52]:


sc.pl.umap(adata_neuronal_disease, color = ["PLP1", ""HLA-DRA", "CXCL9", "cell_death_signature"], cmap="viridis", save="umap_neuronal_gut_atlas_disease_PLP1.png", add_outline=True,
         size=50, alpha=0.7, vmax=3.5)


# In[39]:


# export requirements file
get_ipython().system('pip freeze > requirements_gut_cell_atlas_EGC.txt')


# In[ ]:





# In[ ]:




