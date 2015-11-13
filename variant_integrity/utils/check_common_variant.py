def check_common_variant(genotypes):
    """Check if all individuals share a variant
    
        Args:
            genotypes(iterable): Genotypes
        
        Returns:
            bool: If the individuals share the variant
    """
    
    for genotype in genotypes:
        if not genotype.has_variant:
            return False
    return True
