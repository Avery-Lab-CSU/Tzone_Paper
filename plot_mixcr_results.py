def plot_tcell_clones(extitle='', outfilename='', prefix='./'):
	import pandas as pd
	import numpy as np
	from matplotlib import pyplot as plt
	import glob
	import sys
	print('Starting TCR by counts')
	sys.stdout.flush()
	clist=['C0','C1','C2','C3','C4','C5','C6','C7','C8']
	slist=[]
	sizefactors=pd.read_csv('size_factors.csv', index_col=0)
	fig,(ax1, ax2)=plt.subplots(2,2, figsize=(16,16))
	print('...\tStarting TCRA')
	sys.stdout.flush()
	fl=sorted(glob.glob(prefix+'/*clonotypes.TRA.txt'))
	for i in range(len(fl)):
		holder=pd.read_csv(fl[i], delim_whitespace=True)
		sample_id=fl[i][0:-19]
		if holder.shape[0]==0:
			holder['Sample ID']=sample_id
			ax1[0].bar(x=sample_id,
					   height=0)
		else:
			holder['Sample ID']=sample_id
			holder=holder[['cloneId','cloneCount','Sample ID']]
			blist=[0]
			holder.sort_values('cloneCount', 
							   ascending=False,
							  inplace=True)
			holder['cloneCount']=holder['cloneCount']/sizefactors.loc[sample_id].values
			for i in range(0,len(holder['cloneCount'])-1):
				blist.append(blist[i]+holder.iloc[i,1])
			holder['bottom']=blist
			ll=np.tile(clist, int(np.floor(holder.shape[0]/9)))
			clrs=np.append(ll,clist[0:holder.shape[0]%9])
			holder['color']=clrs
			ax1[0].bar(x=holder['Sample ID'],
					   height=holder['cloneCount'],
					   bottom=holder['bottom'],
					   color=holder['color'])
			   
	ax1[0].tick_params(labelrotation=90)
	ax1[0].set_title(extitle+' TRA Clonotypes')
	ax1[0].set_ylabel('Clone Normalized Read Count')

	####################################################################################TCB
	print('...\tStarting TCRB')
	sys.stdout.flush()
	fl=sorted(glob.glob(prefix+'/*clonotypes.TRB.txt'))
	for i in range(len(fl)):
		holder=pd.read_csv(fl[i], delim_whitespace=True)
		sample_id=fl[i][0:-19]
		if holder.shape[0]==0:
			holder['Sample ID']=sample_id
			ax1[1].bar(x=sample_id,
					   height=0)
		else:
			holder['Sample ID']=sample_id
			holder=holder[['cloneId','cloneCount','Sample ID']]
			blist=[0]
			holder.sort_values('cloneCount', 
							   ascending=False,
							  inplace=True)
			holder['cloneCount']=holder['cloneCount']/sizefactors.loc[sample_id].values
			for i in range(0,len(holder['cloneCount'])-1):
				blist.append(blist[i]+holder.iloc[i,1])
			holder['bottom']=blist
			ll=np.tile(clist, int(np.floor(holder.shape[0]/9)))
			clrs=np.append(ll,clist[0:holder.shape[0]%9])
			holder['color']=clrs
			ax1[1].bar(x=holder['Sample ID'],
					   height=holder['cloneCount'],
					   bottom=holder['bottom'],
					   color=holder['color'])
			   
	ax1[1].tick_params(labelrotation=90)

	ax1[1].set_title(extitle+' TRB Clonotypes')
	ax1[1].set_ylabel('Clone Normalized Read Count')


	####################################################################################TCG
	print('...\tStarting TCRG')
	sys.stdout.flush()
	fl=sorted(glob.glob(prefix+'/*clonotypes.TRG.txt'))
	for i in range(len(fl)):
		holder=pd.read_csv(fl[i], delim_whitespace=True)
		sample_id=fl[i][0:-19]
		if holder.shape[0]==0:
			holder['Sample ID']=sample_id
			ax2[0].bar(x=sample_id,
					   height=0)
		else:
			holder['Sample ID']=sample_id
			holder=holder[['cloneId','cloneCount','Sample ID']]
			blist=[0]
			holder.sort_values('cloneCount', 
							   ascending=False,
							  inplace=True)
			holder['cloneCount']=holder['cloneCount']/sizefactors.loc[sample_id].values
			for i in range(0,len(holder['cloneCount'])-1):
				blist.append(blist[i]+holder.iloc[i,1])
			holder['bottom']=blist
			ll=np.tile(clist, int(np.floor(holder.shape[0]/9)))
			clrs=np.append(ll,clist[0:holder.shape[0]%9])
			holder['color']=clrs
			ax2[0].bar(x=holder['Sample ID'],
					   height=holder['cloneCount'],
					   bottom=holder['bottom'],
					   color=holder['color'])
			   
	ax2[0].tick_params(labelrotation=90)
	ax2[0].set_title(extitle+' TRG Clonotypes')
	ax2[0].set_ylabel('Clone Normalized Read Count')
	ax2[0].set_xlabel('Sample')


	####################################################################################TCD
	print('...\tStarting TCRD')
	sys.stdout.flush()
	fl=sorted(glob.glob(prefix+'/*clonotypes.TRD.txt'))
	for i in range(len(fl)):
		holder=pd.read_csv(fl[i], delim_whitespace=True)
		sample_id=fl[i][0:-19]
		if holder.shape[0]==0:
			holder['Sample ID']=sample_id
			ax2[1].bar(x=sample_id,
					   height=0)
		else:
			holder['Sample ID']=sample_id
			holder=holder[['cloneId','cloneCount','Sample ID']]
			blist=[0]
			holder.sort_values('cloneCount', 
							   ascending=False,
							  inplace=True)
			holder['cloneCount']=holder['cloneCount']/sizefactors.loc[sample_id].values
			for i in range(0,len(holder['cloneCount'])-1):
				blist.append(blist[i]+holder.iloc[i,1])
			holder['bottom']=blist
			ll=np.tile(clist, int(np.floor(holder.shape[0]/9)))
			clrs=np.append(ll,clist[0:holder.shape[0]%9])
			holder['color']=clrs
			ax2[1].bar(x=holder['Sample ID'],
					   height=holder['cloneCount'],
					   bottom=holder['bottom'],
					   color=holder['color'])
	ax2[1].tick_params(labelrotation=90)
	ax2[1].set_title(extitle+' TRD Clonotypes')
	ax2[1].set_ylabel('Clone Normalized Read Count')
	ax2[1].set_xlabel('Sample')
	print('Saving TCR normalized count plot')
	sys.stdout.flush()
	fig.savefig(outfilename+'_Tcell_counts.pdf')

def plot_tcell_clonefractions(extitle='', outfilename='', prefix=''):
	import pandas as pd
	import numpy as np
	from matplotlib import pyplot as plt
	import glob
	import sys
	print('Starting TCR by fraction')
	sys.stdout.flush()
	clist=['C0','C1','C2','C3','C4','C5','C6','C7','C8']
	slist=[]
	fig,(ax1, ax2)=plt.subplots(2,2, figsize=(16,16))
	print('...\tStarting TCRA')
	sys.stdout.flush()
	fl=sorted(glob.glob(prefix+'*clonotypes.TRA.txt'))
	for i in range(len(fl)):
		holder=pd.read_csv(fl[i], delim_whitespace=True)
		sample_id=fl[i][0:-19]
		if holder.shape[0]==0:
			holder['Sample ID']=sample_id
			ax1[0].bar(x=sample_id,
					   height=0)
		else:
			holder['Sample ID']=sample_id
			holder=holder[['cloneId','cloneFraction','Sample ID']]
			blist=[0]
			holder.sort_values('cloneFraction', 
							   ascending=False,
							  inplace=True)
			for i in range(0,len(holder['cloneFraction'])-1):
				blist.append(blist[i]+holder.iloc[i,1])
			holder['bottom']=blist
			ll=np.tile(clist, int(np.floor(holder.shape[0]/9)))
			clrs=np.append(ll,clist[0:holder.shape[0]%9])
			holder['color']=clrs
			ax1[0].bar(x=holder['Sample ID'],
					   height=holder['cloneFraction'],
					   bottom=holder['bottom'],
					   color=holder['color'])
			   
	ax1[0].tick_params(labelrotation=90)
	ax1[0].set_title(extitle+' TRA Clonotypes')
	ax1[0].set_ylabel('Clone Fraction')

	####################################################################################TCB
	print('...\tStarting TCRB')
	sys.stdout.flush()
	fl=sorted(glob.glob(prefix+'/*clonotypes.TRB.txt'))
	for i in range(len(fl)):
		holder=pd.read_csv(fl[i], delim_whitespace=True)
		sample_id=fl[i][0:-19]
		if holder.shape[0]==0:
			holder['Sample ID']=sample_id
			ax1[1].bar(x=sample_id,
					   height=0)
		else:
			holder['Sample ID']=sample_id
			holder=holder[['cloneId','cloneFraction','Sample ID']]
			blist=[0]
			holder.sort_values('cloneFraction', 
							   ascending=False,
							  inplace=True)
			for i in range(0,len(holder['cloneFraction'])-1):
				blist.append(blist[i]+holder.iloc[i,1])
			holder['bottom']=blist
			ll=np.tile(clist, int(np.floor(holder.shape[0]/9)))
			clrs=np.append(ll,clist[0:holder.shape[0]%9])
			holder['color']=clrs
			ax1[1].bar(x=holder['Sample ID'],
					   height=holder['cloneFraction'],
					   bottom=holder['bottom'],
					   color=holder['color'])
			   
	ax1[1].tick_params(labelrotation=90)

	ax1[1].set_title(extitle+' TRB Clonotypes')
	ax1[1].set_ylabel('Clone Fraction')


	####################################################################################TCG
	print('...\tStarting TCRG')
	sys.stdout.flush()
	fl=sorted(glob.glob(prefix+'/*clonotypes.TRG.txt'))
	for i in range(len(fl)):
		holder=pd.read_csv(fl[i], delim_whitespace=True)
		sample_id=fl[i][0:-19]
		if holder.shape[0]==0:
			holder['Sample ID']=sample_id
			ax2[0].bar(x=sample_id,
					   height=0)
		else:
			holder['Sample ID']=sample_id
			holder=holder[['cloneId','cloneFraction','Sample ID']]
			blist=[0]
			holder.sort_values('cloneFraction', 
							   ascending=False,
							  inplace=True)
			for i in range(0,len(holder['cloneFraction'])-1):
				blist.append(blist[i]+holder.iloc[i,1])
			holder['bottom']=blist
			ll=np.tile(clist, int(np.floor(holder.shape[0]/9)))
			clrs=np.append(ll,clist[0:holder.shape[0]%9])
			holder['color']=clrs
			ax2[0].bar(x=holder['Sample ID'],
					   height=holder['cloneFraction'],
					   bottom=holder['bottom'],
					   color=holder['color'])
			   
	ax2[0].tick_params(labelrotation=90)
	ax2[0].set_title(extitle+' TRG Clonotypes')
	ax2[0].set_ylabel('Clone Fraction')
	ax2[0].set_xlabel('Sample')


	####################################################################################TCD
	print('...\tStarting TCRD')
	sys.stdout.flush()
	fl=sorted(glob.glob(prefix+'/*clonotypes.TRD.txt'))
	for i in range(len(fl)):
		holder=pd.read_csv(fl[i], delim_whitespace=True)
		sample_id=fl[i][0:-19]
		if holder.shape[0]==0:
			holder['Sample ID']=sample_id
			ax2[1].bar(x=sample_id,
					   height=0)
		else:
			holder['Sample ID']=sample_id
			holder=holder[['cloneId','cloneFraction','Sample ID']]
			blist=[0]
			holder.sort_values('cloneFraction', 
							   ascending=False,
							  inplace=True)
			for i in range(0,len(holder['cloneFraction'])-1):
				blist.append(blist[i]+holder.iloc[i,1])
			holder['bottom']=blist
			ll=np.tile(clist, int(np.floor(holder.shape[0]/9)))
			clrs=np.append(ll,clist[0:holder.shape[0]%9])
			holder['color']=clrs
			ax2[1].bar(x=holder['Sample ID'],
					   height=holder['cloneFraction'],
					   bottom=holder['bottom'],
					   color=holder['color'])
	ax2[1].tick_params(labelrotation=90)
	ax2[1].set_title(extitle+' TRD Clonotypes')
	ax2[1].set_ylabel('Clone Fraction')
	ax2[1].set_xlabel('Sample')
	print('Saving TCR fraction plot')
	sys.stdout.flush()
	fig.savefig(outfilename+'_Tcell_fraction.pdf')

def plot_bcell_clones(extitle='', outfilename='', prefix='./'):
	import pandas as pd
	import numpy as np
	from matplotlib import pyplot as plt
	import glob
	import sys
	print('Starting IG by counts')
	sys.stdout.flush()
	clist=['C0','C1','C2','C3','C4','C5','C6','C7','C8']
	slist=[]
	sizefactors=pd.read_csv('size_factors.csv', index_col=0)
	fig,(ax1, ax2)=plt.subplots(2,2, figsize=(16,16))
	print('...\tStarting IGH')
	sys.stdout.flush()
	fl=sorted(glob.glob(prefix+'/*clonotypes.TRA.txt'))
	for i in range(len(fl)):
		holder=pd.read_csv(fl[i], delim_whitespace=True)
		sample_id=fl[i][0:-19]
		if holder.shape[0]==0:
			holder['Sample ID']=sample_id
			ax1[0].bar(x=sample_id,
					   height=0)
		else:
			holder['Sample ID']=sample_id
			holder=holder[['cloneId','cloneCount','Sample ID']]
			blist=[0]
			holder.sort_values('cloneCount', 
							   ascending=False,
							  inplace=True)
			holder['cloneCount']=holder['cloneCount']/sizefactors.loc[sample_id].values
			for i in range(0,len(holder['cloneCount'])-1):
				blist.append(blist[i]+holder.iloc[i,1])
			holder['bottom']=blist
			ll=np.tile(clist, int(np.floor(holder.shape[0]/9)))
			clrs=np.append(ll,clist[0:holder.shape[0]%9])
			holder['color']=clrs
			ax1[0].bar(x=holder['Sample ID'],
					   height=holder['cloneCount'],
					   bottom=holder['bottom'],
					   color=holder['color'])
			   
	ax1[0].tick_params(labelrotation=90)
	ax1[0].set_title(extitle+' IGH Clonotypes')
	ax1[0].set_ylabel('Clone Normalized Read Count')

	####################################################################################igk
	print('...\tStarting IGK')
	sys.stdout.flush()
	fl=sorted(glob.glob(prefix+'/*clonotypes.IGK.txt'))
	for i in range(len(fl)):
		holder=pd.read_csv(fl[i], delim_whitespace=True)
		sample_id=fl[i][0:-19]
		if holder.shape[0]==0:
			holder['Sample ID']=sample_id
			ax1[1].bar(x=sample_id,
					   height=0)
		else:
			holder['Sample ID']=sample_id
			holder=holder[['cloneId','cloneCount','Sample ID']]
			blist=[0]
			holder.sort_values('cloneCount', 
							   ascending=False,
							  inplace=True)
			holder['cloneCount']=holder['cloneCount']/sizefactors.loc[sample_id].values
			for i in range(0,len(holder['cloneCount'])-1):
				blist.append(blist[i]+holder.iloc[i,1])
			holder['bottom']=blist
			ll=np.tile(clist, int(np.floor(holder.shape[0]/9)))
			clrs=np.append(ll,clist[0:holder.shape[0]%9])
			holder['color']=clrs
			ax1[1].bar(x=holder['Sample ID'],
					   height=holder['cloneCount'],
					   bottom=holder['bottom'],
					   color=holder['color'])
			   
	ax1[1].tick_params(labelrotation=90)

	ax1[1].set_title(extitle+' IGK Clonotypes')
	ax1[1].set_ylabel('Clone Normalized Read Count')


	####################################################################################IGL
	print('...\tStarting IGL')
	sys.stdout.flush()
	fl=sorted(glob.glob(prefix+'/*clonotypes.IGL.txt'))
	for i in range(len(fl)):
		holder=pd.read_csv(fl[i], delim_whitespace=True)
		sample_id=fl[i][0:-19]
		if holder.shape[0]==0:
			holder['Sample ID']=sample_id
			ax2[0].bar(x=sample_id,
					   height=0)
		else:
			holder['Sample ID']=sample_id
			holder=holder[['cloneId','cloneCount','Sample ID']]
			blist=[0]
			holder.sort_values('cloneCount', 
							   ascending=False,
							  inplace=True)
			holder['cloneCount']=holder['cloneCount']/sizefactors.loc[sample_id].values
			for i in range(0,len(holder['cloneCount'])-1):
				blist.append(blist[i]+holder.iloc[i,1])
			holder['bottom']=blist
			ll=np.tile(clist, int(np.floor(holder.shape[0]/9)))
			clrs=np.append(ll,clist[0:holder.shape[0]%9])
			holder['color']=clrs
			ax2[0].bar(x=holder['Sample ID'],
					   height=holder['cloneCount'],
					   bottom=holder['bottom'],
					   color=holder['color'])
			   
	ax2[0].tick_params(labelrotation=90)
	ax2[0].set_title(extitle+' IGL Clonotypes')
	ax2[0].set_ylabel('Clone Normalized Read Count')
	ax2[0].set_xlabel('Sample')
	fig.delaxes(ax2[1])

	print('Saving IG normalized count plot')
	sys.stdout.flush()
	fig.savefig(outfilename+'_IG_counts.pdf')

def plot_bcell_clonefractions(extitle='', outfilename='', prefix='./'):
	import pandas as pd
	import numpy as np
	from matplotlib import pyplot as plt
	import glob
	import sys

	print('Starting IG by fraction')
	sys.stdout.flush()
	clist=['C0','C1','C2','C3','C4','C5','C6','C7','C8']
	slist=[]
	fig,(ax1, ax2)=plt.subplots(2,2, figsize=(16,16))
	print('...\tStarting IGH')
	sys.stdout.flush()
	fl=sorted(glob.glob(prefix+'/*clonotypes.IGH.txt'))
	for i in range(len(fl)):
		holder=pd.read_csv(fl[i], delim_whitespace=True)
		sample_id=fl[i][0:-19]
		if holder.shape[0]==0:
			holder['Sample ID']=sample_id
			ax1[0].bar(x=sample_id,
					   height=0)
		else:
			holder['Sample ID']=sample_id
			holder=holder[['cloneId','cloneFraction','Sample ID']]
			blist=[0]
			holder.sort_values('cloneFraction', 
							   ascending=False,
							  inplace=True)
			for i in range(0,len(holder['cloneFraction'])-1):
				blist.append(blist[i]+holder.iloc[i,1])
			holder['bottom']=blist
			ll=np.tile(clist, int(np.floor(holder.shape[0]/9)))
			clrs=np.append(ll,clist[0:holder.shape[0]%9])
			holder['color']=clrs
			ax1[0].bar(x=holder['Sample ID'],
					   height=holder['cloneFraction'],
					   bottom=holder['bottom'],
					   color=holder['color'])
			   
	ax1[0].tick_params(labelrotation=90)
	ax1[0].set_title(extitle+' IGH Clonotypes')
	ax1[0].set_ylabel('Clone Fraction')

	####################################################################################IGK
	print('...\tStarting IGK')
	sys.stdout.flush()
	fl=sorted(glob.glob(prefix+'/*clonotypes.IGK.txt'))
	for i in range(len(fl)):
		holder=pd.read_csv(fl[i], delim_whitespace=True)
		sample_id=fl[i][0:-19]
		if holder.shape[0]==0:
			holder['Sample ID']=sample_id
			ax1[1].bar(x=sample_id,
					   height=0)
		else:
			holder['Sample ID']=sample_id
			holder=holder[['cloneId','cloneFraction','Sample ID']]
			blist=[0]
			holder.sort_values('cloneFraction', 
							   ascending=False,
							  inplace=True)
			for i in range(0,len(holder['cloneFraction'])-1):
				blist.append(blist[i]+holder.iloc[i,1])
			holder['bottom']=blist
			ll=np.tile(clist, int(np.floor(holder.shape[0]/9)))
			clrs=np.append(ll,clist[0:holder.shape[0]%9])
			holder['color']=clrs
			ax1[1].bar(x=holder['Sample ID'],
					   height=holder['cloneFraction'],
					   bottom=holder['bottom'],
					   color=holder['color'])
			   
	ax1[1].tick_params(labelrotation=90)

	ax1[1].set_title(extitle+' IGK Clonotypes')
	ax1[1].set_ylabel('Clone Fraction')


	####################################################################################IGL
	print('...\tStarting IGL')
	sys.stdout.flush()
	fl=sorted(glob.glob(prefix+'/*clonotypes.TRG.txt'))
	for i in range(len(fl)):
		holder=pd.read_csv(fl[i], delim_whitespace=True)
		sample_id=fl[i][0:-19]
		if holder.shape[0]==0:
			holder['Sample ID']=sample_id
			ax2[0].bar(x=sample_id,
					   height=0)
		else:
			holder['Sample ID']=sample_id
			holder=holder[['cloneId','cloneFraction','Sample ID']]
			blist=[0]
			holder.sort_values('cloneFraction', 
							   ascending=False,
							  inplace=True)
			for i in range(0,len(holder['cloneFraction'])-1):
				blist.append(blist[i]+holder.iloc[i,1])
			holder['bottom']=blist
			ll=np.tile(clist, int(np.floor(holder.shape[0]/9)))
			clrs=np.append(ll,clist[0:holder.shape[0]%9])
			holder['color']=clrs
			ax2[0].bar(x=holder['Sample ID'],
					   height=holder['cloneFraction'],
					   bottom=holder['bottom'],
					   color=holder['color'])
			   
	ax2[0].tick_params(labelrotation=90)
	ax2[0].set_title(extitle+' IGL Clonotypes')
	ax2[0].set_ylabel('Clone Fraction')
	ax2[0].set_xlabel('Sample')
	fig.delaxes(ax2[1])

	print('Saving IG fraction plot')
	sys.stdout.flush()
	fig.savefig(outfilename+'_IG_fraction.pdf')

if __name__=='__main__':
	import argparse
	parser = argparse.ArgumentParser(description='Plot MIXCR data')
	parser.add_argument('-x', '--prefix', type=str, required=False,help='Path prefix to mixcr output (default: ./)')
	parser.add_argument('-n', '--name', type=str,required=False,help="Name prefix for plots (default: \'\')")
	parser.add_argument('-f', '--file', type=str, required=False,help="Output file name prefix (default: \'\')")
	args = parser.parse_args()
	if args.name and args.prefix and args.file:
		plot_tcell_clones(args.name, args.file, args.prefix)
		plot_tcell_clonefractions(args.name,  args.file, args.prefix)
		plot_bcell_clones(args.name, args.file, args.prefix)
		plot_bcell_clonefractions(args.name, args.file, args.prefix)
	elif args.name and args.prefix:
		plot_tcell_clones(args.name, prefix=args.prefix)
		plot_tcell_clonefractions(args.name, prefix=args.prefix)
		plot_bcell_clones(args.name, prefix=args.prefix)
		plot_bcell_clonefractions(args.name, prefix=args.prefix)
	elif args.name and args.file:
		plot_tcell_clones(args.name, args.file)
		plot_tcell_clonefractions(args.name, args.file)
		plot_bcell_clones(args.name, args.file)
		plot_bcell_clonefractions(args.name, args.file)
	elif args.prefix and args.file:
		plot_tcell_clones(outfilename=args.file, prefix=args.prefix)
		plot_tcell_clonefractions(outfilename=args.file, prefix=args.prefix)
		plot_bcell_clones(outfilename=args.file, prefix=args.prefix)
		plot_bcell_clonefractions(outfilename=args.file, prefix=args.prefix)
	elif args.name:
		plot_tcell_clones(args.name)
		plot_tcell_clonefractions(args.name)
		plot_bcell_clones(args.name)
		plot_bcell_clonefractions(args.name)
	elif args.prefix:
		plot_tcell_clones(prefix=args.prefix)
		plot_tcell_clonefractions(prefix=args.prefix)
		plot_bcell_clones(prefix=args.prefix)
		plot_bcell_clonefractions(prefix=args.prefix)
	elif args.file:
		plot_tcell_clones(outfilename=args.file)
		plot_tcell_clonefractions(outfilename=args.file)
		plot_bcell_clones(outfilename=args.file)
		plot_bcell_clonefractions(outfilename=args.file)
	else:
		plot_tcell_clones()
		plot_tcell_clonefractions()
		plot_bcell_clones()
		plot_bcell_clonefractions()
	