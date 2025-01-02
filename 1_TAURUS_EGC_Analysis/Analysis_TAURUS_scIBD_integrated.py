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
get_ipython().system('wget -o ./vasc_final.h5ad https://zenodo.org/records/14007626/files/vasc_final.h5ad?download=1')

# Load the h5ad file
adata = sc.read_h5ad("./vasc_final.h5ad?download=1")

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


# List all metadata columns
print(adata.obs.columns)


# In[ ]:


# Display rows in .var where 'mt' is True
print(adata.var[adata.var["mt"] == True])

# Display rows in .var where 'rp' is True
print(adata.var[adata.var["rp"] == True])


# In[4]:


# View the first few rows of metadata
print(adata.obs.head())


# In[5]:


# access unique values in a metadata column
display(adata.obs["Ethnicity"].value_counts())
display(adata.obs["Site"].value_counts())
display(adata.obs["Treatment"].value_counts())
display(adata.obs["Remission_status"].value_counts())


# In[6]:


adata.obs["major"].unique()


# In[7]:


# count the number of cells with a certain metadata affiliation
display(adata.obs["major"].value_counts())
display(adata.obs["minor"].value_counts())


# In[8]:


# Compute and plot UMAP
sc.pp.neighbors(adata)
sc.tl.umap(adata)


# In[9]:


# Subset to only cells with metadata value 'Glial' and cells with the umap1 < 2.9 and umap2 < 5
adata_glial = adata[adata.obs["major"] == "Glial"]


# Check the unique values in the 'cell_type' column after subsetting to ensure that subsetting worked
adata_glial.obs["major"].unique()


# In[ ]:


# view basic qc metrics
sc.pl.violin(adata_glial, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_rp'], 
             jitter=0.4, multi_panel=True)


# In[ ]:


# view basic qc metrics
sc.pl.violin(adata_glial, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_rp'], 
             jitter=0.4, multi_panel=True)

upper_lim = np.quantile(adata_glial.obs.n_genes_by_counts.values, .98)


# In[ ]:


upper_lim


# In[ ]:


# filter the cells with high genes by count values
adata_glial = adata_glial[adata_glial.obs.n_genes_by_counts < upper_lim]


# In[ ]:


# filter the cells with high mt percentages
adata_glial = adata_glial[adata_glial.obs.pct_counts_mt < 30]


# In[ ]:


# filter the cells with high ribosomal RNA gene expression
adata_glial = adata_glial[adata_glial.obs.pct_counts_rp < 17.5]


# In[ ]:


# view basic qc metrics
sc.pl.violin(adata_glial, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_rp'], 
             jitter=0.4, multi_panel=True)


# In[10]:


# Create the UMAP plot
sc.pl.umap(
    adata_glial,
    color="major",
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



# In[11]:


# some cells are outliers from the core umap and need to be removed
import pandas as pd

# Extract UMAP coordinates
umap_adata_glia_df = pd.DataFrame(adata_glial.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata_glial.obs.index)

# Add metadata (e.g., cell type, condition)
umap_adata_glia_df["major"] = adata_glial.obs["major"]
umap_adata_glia_df["minor"] = adata_glial.obs["minor"]
umap_adata_glia_df["Disease"] = adata_glial.obs["Disease"]
umap_adata_glia_df["Inflammation"] = adata_glial.obs["Inflammation"]
umap_adata_glia_df["Disease"] = adata_glial.obs["Disease"]




# In[12]:


# to identify the umap coordinates of these cells, use plotly
import plotly.express as px
import plotly.io as pio

# Set the renderer to 'browser'
pio.renderers.default = "browser"

# Create a Plotly scatter plot
fig = px.scatter(
    umap_adata_glia_df,
    x="UMAP1",
    y="UMAP2",
    color="Inflammation",  # Color by inflammation
    hover_data=["Disease"],  # Display condition on hover
    title="Interactive UMAP Plot"
)

fig.show()


# In[13]:


# Check if indices match. If false, reinded the df
print(adata_glial.obs.index.equals(umap_adata_glia_df.index))


# In[14]:


# Reindex the DataFrame to align with adata_glial.obs
umap_adata_glia_df = umap_adata_glia_df.reindex(adata_glial.obs.index)

# Subset the AnnData object
adata_glial = adata_glial[
    (umap_adata_glia_df["UMAP1"] <= 2.9),
    :
].copy()



# In[15]:


# Reindex the DataFrame to align with adata_glial.obs
umap_adata_glia_df = umap_adata_glia_df.reindex(adata_glial.obs.index)

adata_glial = adata_glial[
    (umap_adata_glia_df["UMAP2"] <= 5),
    :
].copy()


# In[16]:


# Reindex the DataFrame to align with adata_glial.obs
umap_adata_glia_df = umap_adata_glia_df.reindex(adata_glial.obs.index)

adata_glial = adata_glial[
    (umap_adata_glia_df["UMAP2"] >= -3.5),
    :
].copy()


# In[47]:


print(f"Number of cells in the subset: {adata_glial.n_obs}")


# In[17]:


# check the umap of the subsetted object
sc.pl.umap(adata_glial, color="major", add_outline=True, legend_loc="on data")



# In[18]:


# Louvain clustering
sc.tl.leiden(adata_glial, resolution=0.1)
sc.pl.umap(adata_glial, color="leiden", legend_loc="on data", size=40, add_outline=True, save="TAURUS_EGC_Clustering.png")  # Visualize clusters


# In[19]:


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
de_table.to_csv("de_per_cluster_EGC_TAURUS_scIBD.csv", index=False)

print("Differential expression results saved to 'de_results_per_cluster.csv'")



# In[21]:


# Generate violin plots for top 20 DE genes
sc.pl.rank_genes_groups_violin(adata_glial, n_genes=20)


# In[22]:


# Heatmap
sc.pl.heatmap(adata_glial, var_names=["PLP1", "S100B", "MPZ", "SPP1", "MAL", "NRXN1", "CD74", "STAT1", "SOCS3", "CXCL9"], 
              groupby=["Inflammation", "Disease"], vmax=8, dendrogram=True)

# Dotplot
sc.pl.dotplot(adata_glial, var_names=["PLP1", "S100B", "SOX10", "MAL", "MPZ", "CCK", "CD74", "DCN", "SPP1", "GFRA3", "NRXN1", "CXCL9"], 
              groupby="leiden")


# In[23]:


# load in cell death list
cell_death = pd.read_csv("cell_death_combined_Hs_scIBD.csv")


# In[25]:


# df to list
cell_death_combined_list = cell_death["Genes"].tolist()


# In[26]:


# Extract expression data from the adata 
# Check which genes from the list are present in adata.var
genes_in_adata = [gene for gene in cell_death_combined_list if gene in adata_glial.var_names]

# Extract expression data for these genes
expression_data = adata_glial[:, genes_in_adata].X.toarray()  # Converts to dense format if sparse


# In[27]:


import numpy as np

# Calculate the mean expression across selected genes for each cell
mean_expression = np.mean(expression_data, axis=1)  # Axis 1: across genes


# In[28]:


# add mean expression as metadata 
adata_glial.obs["cell_death_signature"] = mean_expression




# In[29]:


# List all metadata columns
print(adata_glial.obs.columns)



# In[30]:


# subset to create disease specific objects
adata_EGC_UC = adata_glial[adata_glial.obs["Disease"] == "UC"]

adata_EGC_CD = adata_glial[adata_glial.obs["Disease"] == "CD"]

adata_EGC_CTRL = adata_glial[adata_glial.obs["Disease"] == "Healthy"]



# In[32]:


adata_EGC_UC.obs["Inflammation"].value_counts()


# In[33]:


adata_EGC_CD.obs["Inflammation"].value_counts()


# In[34]:


adata_EGC_CTRL.obs["Inflammation"].value_counts()


# In[39]:


sc.pl.umap(adata_EGC_CTRL, color = ["PLP1", "CXCL9", "cell_death_signature", "Inflammation"], cmap="viridis", save="umap_CTRL_EGCs_markers.png",add_outline=True,
          size=50, alpha=0.7, vmax=3)


# In[40]:


sc.pl.umap(adata_EGC_CTRL, color = ["PLP1", "CXCL9", "cell_death_signature", "Inflammation"], cmap="viridis", save="umap_CTRL_EGCs_markers_PLP.png",add_outline=True,
          size=50, alpha=0.7, vmax=25)


# In[41]:


sc.pl.umap(adata_EGC_CD, color = ["PLP1", "CXCL9", "cell_death_signature", "Inflammation"], cmap="viridis", save="umap_CD_EGCs_markers.png", add_outline=True,
         size=50, alpha=0.7, vmax=3)


# In[42]:


sc.pl.umap(adata_EGC_CD, color = ["PLP1", "CXCL9", "cell_death_signature", "Inflammation"], cmap="viridis", save="umap_CD_EGCs_markers_PLP.png", add_outline=True,
         size=50, alpha=0.7, vmax=25)


# In[43]:


sc.pl.umap(adata_EGC_UC, color = ["PLP1", "CXCL9", "cell_death_signature", "Inflammation"], cmap="viridis", save="umap_UC_EGCs_markers.png", add_outline=True,
         size=50, alpha=0.7, vmax=3)


# In[44]:


sc.pl.umap(adata_EGC_UC, color = ["PLP1", "CXCL9", "cell_death_signature", "Inflammation"], cmap="viridis", save="umap_UC_EGCs_markers_PLP.png", add_outline=True,
         size=50, alpha=0.7, vmax=25)


# In[ ]:


adata_glial.obs["Remission_status"].unique()


# In[ ]:


# split violin plots by remission split by inflammation
sc.pl.violin(adata_glial,
            keys=["PLP1", "CXCL9", "cell_death_signature"],
            groupby="Inflammation",
            jitter=True,
            stripplot=True,
            rotation=45)


# In[1]:


# export requirements file
get_ipython().system('echo "# Python version: $(python --version)" > requirements_TAURUS.txt')
get_ipython().system('pip freeze >> requirements_TAURUS.txt')


# In[ ]:




