## TO DO
# Get headers
# Rewrite bottom code to check if completed
# Set up output dir

import nilearn
from nilearn import datasets
from nilearn import image
from nilearn.maskers import NiftiSpheresMasker
from nilearn.connectome import ConnectivityMeasure
import numpy as np
from nilearn import plotting
import numpy as np
import pandas as pd
#from openpyxl import load_workbook
import matplotlib.pyplot as plt
import sys
import os


def gen_matrices(data, atlas, atlas_networks, atlas_labels, outdir, subj, sesh):


	masker = NiftiSpheresMasker(
		seeds=atlas_labels,
		standardize="zscore",
		memory="nilearn_cache",)

	# Only calculate corr matrix if it doesn't exist
	if not os.path.exists(outdir+subj+"/"+sesh+"/"+subj+"_"+sesh+"_matseitz.csv"):
	##Generate correlation matrices
		masker.fit_transform(data)
		time_series = masker.transform(data) #Still have to add confounds
		correlation_measure = ConnectivityMeasure(kind="correlation")
		correlation_matrix = correlation_measure.fit_transform([time_series])[0]
		np.savetxt(outdir+subj+"/"+sesh+"/"+subj+"_"+sesh+"_matseitz.csv",correlation_matrix, delimiter=",")

	else:
		correlation_matrix = np.loadtxt(outdir+subj+"/"+sesh+"/"+subj+"_"+sesh+"_matseitz.csv", delimiter=",")

	return(masker, correlation_matrix)



def plotting_connectome(correlation_matrix, coords, olddir, subj, sesh):
	np.fill_diagonal(correlation_matrix, 0)
	nilearn.plotting.plot_connectome(
		correlation_matrix, coords, 
		edge_threshold="80%", colorbar=True,
		output_file=olddir+subj+"-"+sesh+"_connectome.png")	

def plotting_matrix(correlation_matrix, atlas_labels, outdir, subj, sesh):

	print(outdir+subj+"/"+subj+"_"+sesh+"_matrixnoconfounds.png")	

	#Fill Diagonal with 0's for visualization
	np.fill_diagonal(correlation_matrix, 0)

	display = nilearn.plotting.plot_matrix(
		correlation_matrix,
		figure=(10,8),
		labels=atlas_labels,
		vmax=0.8,
		vmin=-0.8,
		title="No Confounds",
		reorder=True,
		auto_fit=True
	)

	display.figure.savefig(outdir+subj+"/"+subj+"_"+sesh+"_matrixnoconfounds.png")	

	display.figure.close()

	



def inter_net_conn(masker, correlation_matrix, atlas_networks, network, subj, sesh):

	###Calculate intra-network connectivities
	cormat_df = pd.DataFrame(correlation_matrix,index=atlas_networks,columns=atlas_networks)

	cormat_df__row_Aud = cormat_df.loc[[network],:]
	cormat_df_Aud = cormat_df__row_Aud.loc[:,['Auditory']]
	x = cormat_df_Aud.to_numpy()
	y = np.reshape(x,len(x)*len(x[0]))
	y= y[y!=0]
	mean_cor_Aud = np.mean(y)


	cormat_df__row_Cing = cormat_df.loc[[network],:]
	cormat_df_Cing = cormat_df__row_Cing.loc[:,['CinguloOpercular']]
	#x = np.tril(cormat_df_Cing)
	#y = np.reshape(x,len(x)*len(x))
	x = cormat_df_Cing.to_numpy()
	y = np.reshape(x,len(x)*len(x[0]))
	y= y[y!=0]
	mean_cor_Cing = np.mean(y)

	cormat_df__row_Def = cormat_df.loc[[network],:]
	cormat_df_Def = cormat_df__row_Def.loc[:,['DefaultMode']]
	#x = np.tril(cormat_df_Def)
	#y = np.reshape(x,len(x)*len(x))
	x = cormat_df_Def.to_numpy()
	y = np.reshape(x,len(x)*len(x[0]))
	y= y[y!=0]
	mean_cor_Def = np.mean(y)

	cormat_df__row_DorsAtt = cormat_df.loc[[network],:]
	cormat_df_DorsAtt = cormat_df__row_DorsAtt.loc[:,['DorsalAttention']]
	#x = np.tril(cormat_df_DorsAtt)
	#y = np.reshape(x,len(x)*len(x))
	x = cormat_df_DorsAtt.to_numpy()
	y = np.reshape(x,len(x)*len(x[0]))
	y= y[y!=0]
	mean_cor_DorsAtt = np.mean(y)

	cormat_df__row_FrontPar = cormat_df.loc[[network],:]
	cormat_df_FrontPar = cormat_df__row_FrontPar.loc[:,['FrontoParietal']]
	#x = np.tril(cormat_df_FrontPar)
	#y = np.reshape(x,len(x)*len(x))
	x = cormat_df_FrontPar.to_numpy()
	y = np.reshape(x,len(x)*len(x[0]))
	y= y[y!=0]
	mean_cor_FrontPar = np.mean(y)

	cormat_df__row_MedTemp = cormat_df.loc[[network],:]
	cormat_df_MedTemp = cormat_df__row_MedTemp.loc[:,['MedialTemporalLobe']]
	#x = np.tril(cormat_df_MedTemp)
	#y = np.reshape(x,len(x)*len(x))
	x = cormat_df_MedTemp.to_numpy()
	y = np.reshape(x,len(x)*len(x[0]))
	y= y[y!=0]
	mean_cor_MedTemp = np.mean(y)

	cormat_df__row_ParMed = cormat_df.loc[[network],:]
	cormat_df_ParMed = cormat_df__row_ParMed.loc[:,['ParietoMedial']]
	#x = np.tril(cormat_df_ParMed)
	#y = np.reshape(x,len(x)*len(x))
	x = cormat_df_ParMed.to_numpy()
	y = np.reshape(x,len(x)*len(x[0]))
	y= y[y!=0]
	mean_cor_ParMed = np.mean(y)

	cormat_df__row_Rew = cormat_df.loc[[network],:]
	cormat_df_Rew = cormat_df__row_Rew.loc[:,['Reward']]
	#x = np.tril(cormat_df_Rew)
	#y = np.reshape(x,len(x)*len(x))
	x = cormat_df_Rew.to_numpy()
	y = np.reshape(x,len(x)*len(x[0]))
	y= y[y!=0]
	mean_cor_Rew = np.mean(y)

	cormat_df__row_Sal = cormat_df.loc[[network],:]
	cormat_df_Sal = cormat_df__row_Sal.loc[:,['Salience']]
	#x = np.tril(cormat_df_Sal)
	#y = np.reshape(x,len(x)*len(x))
	x = cormat_df_Sal.to_numpy()
	y = np.reshape(x,len(x)*len(x[0]))
	y= y[y!=0]
	mean_cor_Sal = np.mean(y)

	cormat_df__row_SomDor = cormat_df.loc[[network],:]
	cormat_df_SomDor = cormat_df__row_SomDor.loc[:,['SomatomotorDorsal']]
	#x = np.tril(cormat_df_SomDor)
	#y = np.reshape(x,len(x)*len(x))
	x = cormat_df_SomDor.to_numpy()
	y = np.reshape(x,len(x)*len(x[0]))
	y= y[y!=0]
	mean_cor_SomDor = np.mean(y)

	cormat_df__row_SomLat = cormat_df.loc[[network],:]
	cormat_df_SomLat = cormat_df__row_SomLat.loc[:,['SomatomotorLateral']]
	#x = np.tril(cormat_df_SomLat)
	#y = np.reshape(x,len(x)*len(x))
	x = cormat_df_SomLat.to_numpy()
	y = np.reshape(x,len(x)*len(x[0]))
	y= y[y!=0]
	mean_cor_SomLat = np.mean(y)

	cormat_df__row_Vent = cormat_df.loc[[network],:]
	cormat_df_Vent = cormat_df__row_Vent.loc[:,['VentralAttention']]
	#x = np.tril(cormat_df_Vent)
	#y = np.reshape(x,len(x)*len(x))
	x = cormat_df_Vent.to_numpy()
	y = np.reshape(x,len(x)*len(x[0]))
	y= y[y!=0]
	mean_cor_Vent = np.mean(y)

	cormat_df__row_Vis = cormat_df.loc[[network],:]
	cormat_df_Vis = cormat_df__row_Vis.loc[:,['Visual']]
	#x = np.tril(cormat_df_Vis)
	#y = np.reshape(x,len(x)*len(x))
	x = cormat_df_Vis.to_numpy()
	y = np.reshape(x,len(x)*len(x[0]))
	y= y[y!=0]
	mean_cor_Vis = np.mean(y)

	#Save results to a df
	corvals_ar = np.c_[mean_cor_Aud,mean_cor_Cing,mean_cor_Def,mean_cor_DorsAtt,mean_cor_FrontPar,mean_cor_MedTemp,mean_cor_ParMed,mean_cor_Rew,mean_cor_Sal,mean_cor_SomDor,mean_cor_SomLat,mean_cor_Vent,mean_cor_Vis]
	column = ["mean_cor_Aud","mean_cor_Cing","mean_cor_Def","mean_cor_DorsAtt","mean_cor_FrontPar","mean_cor_MedTemp","mean_cor_ParMed","mean_cor_Rew","mean_cor_Sal","mean_cor_SomDor","mean_cor_SomLat","mean_cor_Vent","mean_cor_Vis"]
	df=pd.DataFrame(corvals_ar, columns=column)
	print(df)

	return df

	


def print_results():

	# Print out the result
	sessions = np.array(sessions)
	#corvals = np.array(FSA_score)
	corvals_ar = np.c_[mean_cor_Aud,mean_cor_Cing,mean_cor_DorsAtt,mean_cor_FrontPar,mean_cor_MedTemp,mean_cor_Vis]
	#np.savetxt('/nethome/jrubio13/FSA/FSA_Results/FSA_results_' +cohort+'_'+phase+'_'+GSR+'.csv', FSA_score_ar, fmt="%s",delimiter=",")
	print(sessions[0]+','+GSR+','+phase+ ","+run+ str(corvals_ar) )



def main():

	atlas = datasets.fetch_coords_seitzman_2018()
	# Loading atlas image stored in 'maps'
	atlas_networks = atlas["networks"]
	# Loading atlas data stored in 'labels'
	atlas_labels = atlas["rois"]
	#print(atlas_labels)


	indir = '/ix1/ginger/dsarpal/lab/reorg/projects/prisma/preproc_bids_2mm/func/'
	outdir = '/ix1/ginger/dsarpal/lab/reorg/projects/prisma/matrix2/'


	for network in atlas_networks:
	# Iter through each matrix
		writer = pd.ExcelWriter("inter_network_{0}.xlsx".format(network), engine='xlsxwriter')

		for subj in os.listdir(indir):
		#define names of input by scan type and run
			if not os.path.exists(outdir+subj):
				os.mkdir(outdir+subj)

			for sesh in os.listdir(indir+subj):
				if not os.path.exists(outdir+subj+"/"+sesh):
					os.mkdir(outdir+subj+"/"+sesh)


				restfile=indir+subj+"/"+sesh+"/brnaswudktm_both_4.nii.gz"

				data = nilearn.image.load_img(restfile)
				masker, correlation_matrix = gen_matrices(data, atlas, atlas_networks, atlas_labels, outdir, subj, sesh)

			
			
				#plotting_matrix(correlation_matrix, atlas_labels, outdir, subj, sesh)
				df = inter_net_conn(masker, correlation_matrix, atlas_networks, network, subj, sesh)
				df.to_excel(writer, sheet_name='{0}_{1}'.format(subj, sesh))

		writer.save()


def run_all_networks(restfile):
	
	finfo = re.search('sub-(?P<subj>[^_/].*ses-(?P<ses>[0-9]{8}',restfile)

	for network in atlas_networks:
	# Iter through each matrix
		writer = pd.ExcelWriter("inter_network_{0}.xlsx".format(network), engine='xlsxwriter')

		for subj in os.listdir(indir):
		#define names of input by scan type and run
			if not os.path.exists(outdir+subj):
				os.mkdir(outdir+subj)

			for sesh in os.listdir(indir+subj):
				if not os.path.exists(outdir+subj+"/"+sesh):
					os.mkdir(outdir+subj+"/"+sesh)


				restfile=indir+subj+"/"+sesh+"/brnaswudktm_both_4.nii.gz"

				data = nilearn.image.load_img(restfile)
				masker, correlation_matrix = gen_matrices(data, atlas, atlas_networks, atlas_labels, outdir, subj, sesh)

			


atlas = datasets.fetch_coords_seitzman_2018()
# Loading atlas image stored in 'maps'
atlas_networks = atlas["networks"]
# Loading atlas data stored in 'labels'
atlas_labels = atlas["rois"]
#print(atlas_labels)

preproc_files=glob.glob("/path/to/sub*/*/brnaswudktm_both_4.nii.gz)


main()
