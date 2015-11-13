def check_high_quality(genotypes, gq_treshold):
    """Check if all genotype calls are above treshold
    
        Args:
            genotypes (iterable): an iterator with Genotypes
            gq_treshold (int): The genotype treshold
        
        Returns:
            bool: If gq treshold was met or not
    """
    
    for genotype in genotypes:
        if genotype.genotype_quality < gq_treshold:
            return False
    
    return True
    
