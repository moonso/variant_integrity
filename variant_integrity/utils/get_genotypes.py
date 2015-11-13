from vcftoolbox import Genotype

def get_genotypes(variant, analysis_individuals):
    """
    Return Genotype objects for the analysis individuals
    
        Args:
            variant (dict): A variant dictionary
            analysis_individuals (iterator(str)): The individual ids 
        
        Returns:
            genotypes (dict): A dictionary with individual ids as keys and 
                                Genotype objects as values
    """
    genotypes = {}
    gt_format = variant['FORMAT'].split(':')
    for ind_id in analysis_individuals:
        genotypes[ind_id] = Genotype(**dict(zip(
                    gt_format,
                    variant[ind_id].split(':')
                    )))
    return genotypes
    