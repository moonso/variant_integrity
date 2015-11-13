def check_mendelian_error(child, mother, father):
    """Check if there is a mendelian error for a trio
    
        Args:
            child (Genotype)
            mother (Genotype)
            father (Genotype)
        
        Returns:
         bool : if it is a mendelian error or not 
    """
    mendelian_error = False
    
    # If child is homozygote alternative both parents has to have the variant
    if child.homo_alt:
        if not (mother.has_variant and father.has_variant):
            return True
    else:
        if not (mother.has_variant or father.has_variant):
            return True
    
    return False
    