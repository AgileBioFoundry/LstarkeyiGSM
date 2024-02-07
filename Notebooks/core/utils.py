from matplotlib import pyplot as plt
from matplotlib import colors
import numpy as np
import pandas as pd
import cobra
import seaborn as sns
from sklearn.metrics import confusion_matrix, matthews_corrcoef


def background_gradient(s, cmap='seismic', text_color_threshold=0.408):
	'''
	obtain the min and max values of growth and metabolic function from experimental data. 
	
	Args: 
		s
		
		camp:
		
		
		text_color_threshold
	'''
	# obtain the min/max of the data.
	lim = max(abs(s.min().min()),abs(s.max().max()))
	
	rng = 2.0*lim 
	
	norm = colors.Normalize(-lim - (rng * 0.2), lim + (rng * 0.2))
	rgbas = plt.cm.get_cmap(cmap)(norm(s.values))
	
	# 
	def relative_luminance(rgba):
		r, g, b = (x / 12.92 if x <= 0.03928 else ((x + 0.055) / 1.055 ** 2.4) for x in rgba[:3])
		return 0.2126 * r + 0.7152 * g + 0.0722 * b
	
	#
	def css(rgba):
		dark = relative_luminance(rgba) < text_color_threshold
		text_color = '#f1f1f1' if dark else '#000000'
		return 'background-color: {b};color: {c};'.format(b=colors.rgb2hex(rgba), c=text_color)

	# 
	if s.ndim == 1:
		return [css(rgba) for rgba in rgbas]
	else:
		return pd.DataFrame([[css(rgba) for rgba in row] for row in rgbas], index=s.index, columns=s.columns)

#	  
def Show_Data(x):
	'''
	plot the data color maps.
	
	Args: 
		x
	'''
	
	display(Fitness.loc[temp].style.apply(background_gradient, cmap='seismic', axis=None).highlight_null('lightgrey'))
	return;

def PerformFBAPredictions(model,Biolog,Biolog_Prediction,Biolog_Media,Biolog_in_model,threshold = 1e-6): 
	with model:

		# set uptake reactions of other carbon and nitrogen sources to zero and then set the biolog assay well source for uptake. 
		model.reactions.get_by_id('ATPM').lower_bound = 0.0
		model.reactions.get_by_id('EX_glc__D_e').lower_bound = 0.0
		model.reactions.get_by_id('EX_nh4_e').lower_bound = 0.0
		model.reactions.get_by_id('EX_pi_e').lower_bound = 0.0
		model.reactions.get_by_id('EX_so4_e').lower_bound = 0.0

		# keep track of the condition being examined.
		for i, row in Biolog.iterrows():
			Biolog_Prediction.loc[i,'PlateType'] = row['PlateType']
			Biolog_Prediction.loc[i,'Experiment'] = row['Experiment']
			Biolog_Prediction.loc[i,'Row'] = row['Row']
			Biolog_Prediction.loc[i,'Column'] = row['Column']
			# Biolog_Prediction.loc[i,'Data'] = row['Average']
			# Biolog_Prediction.loc[i,'590_avg'] = row['AverageValue_590']
			# Biolog_Prediction.loc[i,'750_avg'] = row['AverageValue_750']
			# Biolog_Prediction.loc[i,'750_pval'] = row['p-values_750']
			# Biolog_Prediction.loc[i,'590_pval'] = row['p-values_590']
			Biolog_Prediction.loc[i,'Data_TF'] = row['SigGrowth?']

			# object to keep track of whether metabolite/reaction of source is actually in the model.
			temp = 1

			# check that there are external metabolites and set them to allow for uptake.
			if i in Biolog_Media.index and Biolog_in_model.loc[i,'External']:
				for x in Biolog_Media.loc[i][1:]:

					# Check if there are multiple nutrient sources in the well, set both to uptake.
					if ',' in x:
						for y in x.split(','):

							# statement to double check that reactions are in model. 
							try:
								model.reactions.get_by_id(y).lower_bound = -10.0

							# case where reaction is not in model.
							except:
								print('warning - reaction {} is not in model'.format(x))
								Biolog_Prediction.loc[i,'Prediction_TF'] =np.nan
								Biolog_Prediction.loc[i,'Prediction'] = np.nan
								temp = 0

					# if there is only one nutrient source, set the uptake to -10.
					else:

						# statement to double check that reactions are in model. 
						try:
							model.reactions.get_by_id(x).lower_bound = -10.0

						# case where reaction is not in model.
						except:
							print('warning - reaction {} is not in model'.format(x))
							Biolog_Prediction.loc[i,'Prediction_TF'] =np.nan
							Biolog_Prediction.loc[i,'Prediction'] = np.nan
							temp = 0


				# optimize FBA object.
				sol = model.optimize()			  


				# reset GSM to prevent uptake so we can continue to the next biolog condition.
				for x in Biolog_Media.loc[i][1:]:

					# multiple nutrients being uptaken. 
					if ',' in x:

						# statement to double check that reactions are in model. 
						try:
							for y in x.split(','):
								model.reactions.get_by_id(y).lower_bound = 0.0

						# case when reactions are not in model. 
						except:
							pass

					# cases where only one nutrient is being uptaken.
					else:

						# statement to double check that reactions are in model. 
						try:
							model.reactions.get_by_id(x).lower_bound = 0.0

						# case when reactions are not in model. 
						except:
							pass


				# check if the solution had feasible growth.
				if sol.status == 'optimal':
					# record the predicted growth.
					Biolog_Prediction.loc[i,'Prediction'] = sol.objective_value

					# classify the growth as successful if greater threshold
					Biolog_Prediction.loc[i,'Prediction_TF'] = sol.objective_value > threshold
					i = 0

				# provide warning that GSM did predict growth under biolog condition.
				else:
					print(i, 'non-optimal')
					Biolog_Prediction.loc[i,'Prediction'] = np.nan
					Biolog_Prediction.loc[i,'Prediction_TF'] = np.nan

			# if there are no external metabolites, check if there are internal ones. 
			# create a new export reaction and external metabolite for the intended source.
			elif i in Biolog_Media.index and Biolog_in_model.loc[i,'Internal']:
				continue

				# obtain the reaction name.
				for x in Biolog_Media.loc[i,Biolog_Media.loc[i,'Source']].split(','):

					# check if it is in the model reactions. 
					if not x in model.reactions:

						# create a new reaction using hydrogen exchange reaction as a scaffold. 
						r = model.reactions.get_by_id('EX_h_e').copy()
						r.id = x
						model.add_reactions([r])
						r.add_metabolites({'h_e': 1.0, x.replace('EX_','').rsplit('_',1)[0]+'_c': -1.0})

				# after creating reaction - set sources to allow for growth - same as instances above. 
				for x in Biolog_Media.loc[i][1:]:

					# Check if there are multiple nutrient sources in the well, set both to uptake.
					if ',' in x:
						for y in x.split(','):

							# statement to double check that reactions are in model. 
							try:
								model.reactions.get_by_id(y).lower_bound = -10.0

							# case where reaction is not in model.
							except:
								print('warning - reaction {} is not in model'.format(x))
								Biolog_Prediction.loc[i,'Prediction_TF'] = np.nan
								Biolog_Prediction.loc[i,'Prediction'] = np.nan
								temp = 0

					# if there is only one nutrient source, set the uptake to -10.
					else:

						# statement to double check that reactions are in model. 
						try:
							model.reactions.get_by_id(x).lower_bound = -10.0

						# case where reaction is not in model.
						except:
							print('warning - reaction {} is not in model'.format(x))
							Biolog_Prediction.loc[i,'Prediction_TF'] = np.nan
							Biolog_Prediction.loc[i,'Prediction'] = np.nan
							temp = 0


				# optimize FBA object.
				sol = model.optimize()			  


				# reset GSM to prevent uptake so we can continue to the next biolog condition.
				for x in Biolog_Media.loc[i][1:]:

					# multiple nutrients being uptaken. 
					if ',' in x:

						# statement to double check that reactions are in model. 
						try:
							for y in x.split(','):
								model.reactions.get_by_id(y).lower_bound = 0.0

						# case when reactions are not in model. 
						except:
							pass

					# cases where only one nutrient is being uptaken.
					else:

						# statement to double check that reactions are in model. 
						try:
							model.reactions.get_by_id(x).lower_bound = 0.0

						# case when reactions are not in model. 
						except:
							pass


				# check if the solution had feasible growth.
				if sol.status == 'optimal':
					# record the predicted growth.
					Biolog_Prediction.loc[i,'Prediction'] = sol.objective_value

					# classify the growth as successful if greater threshold
					Biolog_Prediction.loc[i,'Prediction_TF'] = sol.objective_value > threshold
					i = 0

				# provide warning that GSM did predict growth under biolog condition.
				else:
					print(i, 'non-optimal')
					Biolog_Prediction.loc[i,'Prediction'] = np.nan
					Biolog_Prediction.loc[i,'Prediction_TF'] = np.nan
			else:
				Biolog_Prediction.loc[i,'Prediction'] = np.nan
				Biolog_Prediction.loc[i,'Prediction_TF'] = np.nan
				#Biolog_Prediction.loc[i,'Prediction'] = 0.0
				#Biolog_Prediction.loc[i,'Prediction_TF'] = False
	return(Biolog_Prediction)



def gapfil_reactions(model,eco,Biolog_Prediction,Nutrients,Biolog_in_model,threshold=1e-6):

	'''
	This function returns performs FBA gapfilling (https://cobrapy.readthedocs.io/en/latest/gapfilling.html)
	and returns a list of reactions that enable growth on a substrate.
	
	inputs:
		model - fba.model:
			target model to gapfill.
		eco - fba.model
			model that is used as template to gapfill.
		Biolog_Prediction - pd.DataFrame():
			input data that contains the true/false for metabolite growth in model.
		Nutrients - pd.DataFrame():
			a dataframe formation of the type of source for growth.
		Biolog_in_model - pd.DataFrame():
			a dataframe formation of a dict that contains the reaction names in the model.
		threshold - np.float:
			minimal value needed to be achieved to assert successful growth.
			 
	outputs:
		model - fba.model:
			target model that has been gapfilled. 
		rxns_neededToFix_noGrowth - pd.DataFrame():
			reactions needed to be to enable growth. 
	
	'''

	new_reactions = []
	new_reactions_name = []
	new_reactions_gpr = []
	new_reaction_stoich	 = []
	sources_fixed = []
	for nutrient in Nutrients:

		# print the nutrient source being tested.
		print(nutrient,'\n\n')
		
		
		
		# determine biolog conditions that were able to be tested in the model.
		temp = Biolog_Prediction.index[~np.isnan(Biolog_Prediction.Prediction)]

		# experimental data.
		y_data = Biolog_Prediction.Data_TF[temp].astype(int)

		# model prediction.
		y_pred = Biolog_Prediction.Prediction_TF[temp].astype(int)
		
		# get index where the prediction was incorrect in predicting no growth.
		t = Biolog_Prediction.loc[Biolog_Prediction.loc[y_pred[(y_data!=y_pred) & (y_pred==0)].index,:].index,'Experiment']==nutrient
		t = t[t==True]

		# obtain the nutrient failed sources that failed.
		Sources = Biolog_in_model.loc[t.index,:].Exchange.values

		# obtain the ATPM value of the template model.
		ATPM_template = eco.reactions.ATPM.lower_bound

		# iterate through the failed source reactions.
		for source in Sources:

			# set bounds to grow on source. 
			try:

				# allow uptake. 
				model.reactions.get_by_id(source).lower_bound = -10

				# set bounds of default source to zero if used as sole nutrient in model.
				if nutrient == 'Carbon':
					model.reactions.get_by_id('EX_glc__D_e').lower_bound = 0
				elif nutrient == 'Nitrogen':
					model.reactions.get_by_id('EX_nh4_e').lower_bound = 0
				elif nutrient == 'Phosphorus':
					model.reactions.get_by_id('EX_pi_e').lower_bound = 0
				elif nutrient == 'Sulfur':
					model.reactions.get_by_id('EX_so4_e').lower_bound = 0


	#			  print(model.optimize().objective_value, '\t', model)

			except:
				print('warning \t\t', model, 'does not have reaction')

			# set bounds in e. coli model to grow only on source.
			try:

				# allow source uptake.
				eco.reactions.get_by_id(source).lower_bound = -10
			
				# set bounds of default source to zero if used as sole nutrient in model.
				if nutrient == 'Carbon':
					eco.reactions.get_by_id('EX_glc__D_e').lower_bound = 0
				elif nutrient == 'Nitrogen':
					eco.reactions.get_by_id('EX_nh4_e').lower_bound = 0
				elif nutrient == 'Phosphorus':
					eco.reactions.get_by_id('EX_pi_e').lower_bound = 0
				elif nutrient == 'Sulfur':
					eco.reactions.get_by_id('EX_so4_e').lower_bound = 0
	#			  print(eco.optimize().objective_value, '\t', eco)

				print('----\n',source)
					
				eco.reactions.ATPM.lower_bound = 0
				try:
					template_obj = eco.optimize().objective_value
				except:
					print('template model infeasible')
			
				if template_obj>threshold:
			
					# perform gap filling. 
					try:
						gapfiller = cobra.flux_analysis.gapfilling.GapFiller(model, eco, demand_reactions=False, integer_threshold=1e-9)
						gapfiller.model.solver.configuration.tolerances.feasibility = 1e-9
						gapfiller.model.solver.configuration.tolerances.integrality = 1e-9
						# gapfiller.model.solver.configuration.tolerances.optimality = 1e-9
	
						print('gapfiller initialized')
						# save the reactions neeeded
						reacs = gapfiller.fill()
		#				  print('***')

		#				  print(source)
						print('gapfiller succeeded')
						# add reactions to the model.
						for r in reacs[0]:
							new_reactions.append(r.id)
							new_reactions_name.append(r.name)
							new_reactions_gpr.append(r.gene_reaction_rule)
							new_reaction_stoich.append(r.reaction)
							sources_fixed.append(source)
							print(r, '\t', r.name)
							
							# potential to add reactions to the model.
							model.add_reactions([r.copy()])
						print('\n')
				
					except Exception as e:
						print(e)
						try:
							gapfiller.model.solver.configuration.tolerances.feasibility = 1e-5
							gapfiller.model.solver.configuration.tolerances.integrality = 1e-5
	
							# save the reactions neeeded
							reacs = gapfiller.fill()

							print('gapfiller succeeded round 2')
							# add reactions to the model.
							for r in reacs[0]:
								new_reactions.append(r.id)
								new_reactions_name.append(r.name)
								new_reactions_gpr.append(r.gene_reaction_rule)
								new_reaction_stoich.append(r.reaction)
								print(r, '\t', r.name)
								sources_fixed.append(source)
								
								# potential to add reaction to model. 
# 								model.add_reactions([r.copy()])
						except:
							print('gapfiller failed second time')
				else:
					print('template model growth failed ({})'.format(template_obj))
	
			except:
				print('warning \t\t', eco, 'does not have reaction',source)




			#print(reacs)


			# reset all of the sources to no uptake for next iteration.
			try:
				model.reactions.get_by_id(source).lower_bound = 0

				# reset source
				if nutrient == 'Carbon':
					model.reactions.get_by_id('EX_glc__D_e').lower_bound = -10
				elif nutrient == 'Nitrogen':
					model.reactions.get_by_id('EX_nh4_e').lower_bound = -1000
				elif nutrient == 'Phosphorus':
					model.reactions.get_by_id('EX_pi_e').lower_bound = -1000
				elif nutrient == 'Sulfur':
					model.reactions.get_by_id('EX_so4_e').lower_bound = -1000


		#		  model.reactions.get_by_id('EX_glc__D_e').lower_bound = 0
		#		  print(model.optimize().objective_value, '\t', model)
			except:
				print('warning \t\t', model, 'does not have reaction')

			# reset template based model. 
			try:
				eco.reactions.get_by_id(source).lower_bound = 0

				# reset source
				if nutrient == 'Carbon':
					eco.reactions.get_by_id('EX_glc__D_e').lower_bound = -10
				elif nutrient == 'Nitrogen':
					eco.reactions.get_by_id('EX_nh4_e').lower_bound = -1000
				elif nutrient == 'Phosphorus':
					eco.reactions.get_by_id('EX_pi_e').lower_bound = -1000
				elif nutrient == 'Sulfur':
					eco.reactions.get_by_id('EX_so4_e').lower_bound = -1000
				
				eco.reactions.ATPM.lower_bound = ATPM_template

		#		  eco.reactions.get_by_id('EX_glc__D_e').lower_bound = 0
		#		  print(eco.optimize().objective_value, '\t', eco)
			except:
				print('warning \t\t', eco, '- template - does not have reaction')
	rxns_neededToFix_noGrowth = pd.DataFrame([new_reactions, new_reactions_gpr,new_reaction_stoich, new_reactions_name,sources_fixed])
	return(model,rxns_neededToFix_noGrowth)