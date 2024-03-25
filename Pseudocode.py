#Aceder ao obs do ficheiro "adata_final_{d}_cca_features.h5ad" (output do leiden) d= Astrocytes ou Immune

#Após definir qual a melhor resolução para "seriação primária" distinção robusta de clusters "isolados" (Astrocytes == 0.4? / Immune == 1.2?), criar uma coluna nova no obs que servirá para assinalar (a partir dos scores do leiden armazenados no obs) quais as células que pertencem a esses clusters nessa resolução e quais os clusters a que pertencem (tomando o cuidado de alterar manualmente o nome de maneira a que se saiba que pertencem a um cluster da respetiva resolução ex= 04C3). 

#Selecionar as que células que ficam assinaladas como não pertencentes a nenhum desses clusters, e atribuir o valor do leiden score que lhes é atribuido na resolução da "seriação secundária" (Astrocytes == 0.7? / Immune == 0.9?)  


##Fazer um output de umaps da resolução test 1 e 2 

### Execute TODO 1

### Select the cells
# The selected cells MUST belong to BOTH good clusters list (&&)
    # mask = obs[rXX].isin([lXX]) && obs[rYY].isin([lYY])
    # adata.subsetting[mask] 
# The remainig cells are all the others
    # adata.subsetting[~mask] (! or ~, depending on the language)
    # TODO 2: Verify that mask negation works as intendend
# Get cell tag names
    # good_cell_names = adata.obs_names.to_list()   
    # repeat for bad_cells      



### Copy info to new fusion cluster
## Add good_cell info
# adata_cca_feature[good_cell_names, :].obs[fusion_cluster] = .obs['test_leiden_n15_r2']
## Convert cluster names
#for key, value in dict_good_cells.items():
   #adata[   adata.obs[fusion_cluster].isin([key])  ].obs[fusion_cluster] = value
# ...
## Add bad_cell info
# ...
## Convert cluster names
# ....
## Convert fusion cluster to categories:
# ...
# TODO 3: print of categories names
#print(small_adata.obs['fusion_cluster'].cat.categories)

### Output
# Fazer um output de umaps da resolução test_leiden_r1 , r2, r1_2 (fusion of r1 and r2) 

### TODO 1:
## Print cluster elements
# print(adata.obs['test_leiden_n15_r1'])
# print(adata.obs['test_leiden_n15_r1'].cat.categories))
## How to select cells from a cluster?
# Before selsction:
# print(adata)
# small_adata = adata[adata.obs['test_leiden_n15_r1' == '0']]
# print(small_adata)
# print(small_adata.obs['test_leiden_n15_r1'].cat.categories) 
# med_adata = adata[adata.obs['test_leiden_n15_r1'].isin(['0', '1', '10', '15', '24', '30'])]
# print(med_adata)
# print(med_adata.obs['test_leiden_n15_r1'].cat.categories)
### END TODO 1

### How to convert to categories:
# data[d].obs[f'{es}_rank'] = None
# code... code ...
#   data[d].obs.loc[top_cells, f'{es}_rank'] = 'top'
#   data[d].obs.loc[low_cells, f'{es}_rank'] = 'low'
# 
# When column is finished: 
# data[d].obs[f'{es}_rank'] = data[d].obs[f'{es}_rank'].astype(dtype='category')
#adata.obs['fusion_cluster'] = adata.obs['fusion_cluster'].astype("category")












#combined_mask = None
          #for r in Good_clusters_r2:
          #    for i in Good_clusters_r1
          #    mask = adata.obs["test_leiden_n15_r2.00"].isin([r]) && adata.obs["test_leiden_n15_r2.00"]
          #   if combined_mask is None:
          #      combined_mask = mask
          #   else:
          #      combined_mask |= mask
          
          #selected_cells = adata[combined_mask]