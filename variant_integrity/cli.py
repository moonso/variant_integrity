# -*- coding: utf-8 -*-
from __future__ import division

import sys
import os
import logging
import itertools
import json

import click

import variant_integrity

from ped_parser import FamilyParser
from vcftoolbox import (HeaderParser, get_variant_dict, get_variant_id)

from .log import configure_stream, LEVELS
from variant_integrity.utils import (check_individuals, get_genotypes,
check_mendelian_error, check_high_quality, check_common_variant)

logger = logging.getLogger(__name__)


@click.group()
@click.argument('variant_file', 
                type=click.File('r'),
                metavar='<vcf_file> or -'
)
@click.option('-f', '--family_file',
                type=click.File('r'),
                metavar='<ped_file>'
)
@click.option('-t' ,'--family_type', 
                type=click.Choice(['ped', 'alt', 'cmms', 'mip']), 
                default='ped',
                help='If the analysis use one of the known setups, please specify which one.'
)
@click.option('-g' ,'--gq_treshold', 
                default=20,
                help='The genotype quality to consider a variant'
)
@click.option('--to_json',
                is_flag=True,
                help='If json output.'
)
@click.option('-v', '--verbose', 
                count=True, 
                default=2
)
@click.option('-o', '--outfile',
                type=click.File('w'),
                help='Specify the path to a file where results should be stored.'
)
@click.version_option(variant_integrity.__version__)
@click.pass_context
def cli(ctx, variant_file, family_file, family_type, gq_treshold, to_json,
    outfile, verbose):
    """Check for pedigree inconsistensies."""
    # configure root logger to print to STDERR
    loglevel = LEVELS.get(min(verbose, 3))
    configure_stream(level=loglevel)

    if not family_file:
        logger.error("Please provide a family file with -f/--family_file")
        logger.info("Exiting")
        sys.exit(1)
    
    logger.info("Setting up a family parser")
    family_parser = FamilyParser(family_file, family_type)
    logger.debug("Family parser done")
    # The individuals in the ped file must be present in the variant file:
    families = family_parser.families
    logger.info("Families used in analysis: {0}".format(
                    ','.join(list(families.keys()))))
    
    ctx.gq_treshold = gq_treshold
    ctx.to_json = to_json
    ctx.outfile = outfile

    ctx.families = families
    ctx.individuals = family_parser.individuals
    
    head = HeaderParser()
    
    for line in variant_file:
        line = line.rstrip()

        if line.startswith('#'):
            if line.startswith('##'):
                head.parse_meta_data(line)
            else:
                head.parse_header_line(line)
        else:
            break
    #Add the first variant to the iterator
    if line:
        variant_file = itertools.chain([line], variant_file)
    
    try:
        check_individuals(family_parser.individuals, head.individuals)
    except IOError as e:
        logger.error(e)
        logger.info("Individuals in PED file: {0}".format(
                        ', '.join(family_parser.individuals)))
        logger.info("Individuals in VCF file: {0}".format(', '.join(vcf_individuals)))
        logger.info("Exiting...")
        ctx.abort()
    
    ctx.variant_file = variant_file
    ctx.header_line = head.header
    

@cli.command()
@click.version_option(variant_integrity.__version__)
@click.pass_context
def mendel(ctx):
    """Check mendelian errors in all trios"""
    logger.info("Running variant_integrity mendel {0}".format(
        variant_integrity.__version__))
    
    print_columns = ['ind_id', 'fraction_of_errors', 'mendelian_errors', 'number_calls']
    # Children is a dictionary of children that counts the number of errors
    trios = []
    children = {}
    analysis_individuals = set()
    
    for family in ctx.parent.families:
        family_object = ctx.parent.families[family]
        for trio in family_object.trios:
            trio_individuals = {
                'mother':None,
                'father':None,
                'child':None,
            }
            for ind_id in trio:
                analysis_individuals.add(ind_id)
                individual_object = ctx.parent.individuals[ind_id]
                if individual_object.mother in trio:
                    trio_individuals['child'] = ind_id
                elif individual_object.sex == 1:
                    trio_individuals['father'] = ind_id
                else:
                    trio_individuals['mother'] = ind_id
            trios.append(trio_individuals)
            logger.info("Trio found: {0}".format(
                ', '.join(list(trio_individuals.values()))
            ))
    
    logger.info("Individuals included in analysis: {0}".format(
                    ','.join(list(analysis_individuals))))
    
    for trio in trios:
        children[trio['child']] = dict(zip(
            print_columns, [trio['child'], 0, 0, 0]))
    
    
    for line in ctx.parent.variant_file:
        variant_dict = get_variant_dict(
            variant_line=line, 
            header_line=ctx.parent.header_line
        )
        
        logger.debug("Checking genotype calls for variant {0}".format(
            get_variant_id(variant_dict=variant_dict)
        ))
        
        genotypes = get_genotypes(variant_dict, analysis_individuals)
        
        for trio in trios:
            child_id = trio['child']
            child_genotype = genotypes[child_id]
            mother_genotype = genotypes[trio['mother']]
            father_genotype = genotypes[trio['father']]
            trio_genotypes = [
                child_genotype,
                mother_genotype,
                father_genotype
            ]
            #First check if the child has the variant:
            if child_genotype.has_variant:
                # If all individuals are high quality we count the variant
                if check_high_quality(trio_genotypes, ctx.parent.gq_treshold):
                    children[child_id]['number_calls'] += 1
                    if check_mendelian_error(child_genotype, mother_genotype, father_genotype):
                        children[child_id]['mendelian_errors'] += 1
    
    results = []
    
    for child_id in children:
        child_info = children[child_id]
        errors = child_info['mendelian_errors']
        variants = child_info['number_calls']
        percentage = errors/variants
        
        child_info['fraction_of_errors'] = round(percentage, 3)
        results.append(child_info)
    
    to_json = ctx.parent.to_json
    outfile = ctx.parent.outfile
    
    if to_json:
        if outfile:
            json.dump(results, outfile)
        else:
            print(json.dumps(results))
    else:
        if outfile:
            outfile.write("#{0}\n".format('\t'.join(print_columns)))
        else:
            print("#{0}".format('\t'.join(print_columns)))
        
        for result in results:
            print_line = "{0}\t{1}\t{2}\t{3}".format(
                    result['ind_id'], result['fraction_of_errors'],
                    result['mendelian_errors'], result['number_calls']
                )
            if outfile:
                outfile.write("{0}\n".format(print_line))
            else:
                print(print_line)

@cli.command()
@click.version_option(variant_integrity.__version__)
@click.pass_context
def father(ctx):
    """Check number of variants in common"""
    logger.info("Running variant_integrity father version {0}".format(
        variant_integrity.__version__))
    
    print_columns = ['ind_id', 'fraction_of_common_variants', 'common_variants', 'number_calls']
    # Children is a dictionary of children that counts the number of errors
    duos = []
    children = {}
    analysis_individuals = set()
    
    for ind_id in ctx.parent.individuals:
        individual_object = ctx.parent.individuals[ind_id]
        if individual_object.father != '0':
            duo = {
                'child': ind_id,
                'father': individual_object.father
            } 
            
            analysis_individuals.add(ind_id)
            analysis_individuals.add(individual_object.father)
            
            duos.append(duo)
            logger.info("Duo found: {0}".format(
                ', '.join(list(duo.values()))
            ))
    
    logger.info("Individuals included in analysis: {0}".format(
                    ','.join(list(analysis_individuals))))
    
    for duo in duos:
        children[duo['child']] = dict(zip(
            print_columns, [duo['child'], 0, 0, 0]))
    
    
    for line in ctx.parent.variant_file:
        variant_dict = get_variant_dict(
            variant_line=line, 
            header_line=ctx.parent.header_line
        )
        
        logger.debug("Checking genotype calls for variant {0}".format(
            get_variant_id(variant_dict=variant_dict)
        ))
        
        genotypes = get_genotypes(variant_dict, analysis_individuals)
        
        for duo in duos:
            child_id = duo['child']
            child_genotype = genotypes[child_id]
            father_genotype = genotypes[duo['father']]
            duo_genotypes = [
                child_genotype,
                father_genotype
            ]
            #First check if the child has the variant:
            if child_genotype.has_variant:
                # If child have high quality we count the variant
                if check_high_quality([child_genotype], ctx.parent.gq_treshold):
                    children[child_id]['number_calls'] += 1
                    if check_common_variant(duo_genotypes):
                        children[child_id]['common_variants'] += 1
    
    results = []
    
    for child_id in children:
        child_info = children[child_id]
        common = child_info['common_variants']
        variants = child_info['number_calls']
        percentage = common/variants
        
        child_info['fraction_of_common_variants'] = round(percentage, 3)
        results.append(child_info)
    
    to_json = ctx.parent.to_json
    outfile = ctx.parent.outfile
    
    if to_json:
        if outfile:
            json.dump(results, outfile)
        else:
            print(json.dumps(results))
    else:
        if outfile:
            outfile.write("#{0}\n".format('\t'.join(print_columns)))
        else:
            print("#{0}".format('\t'.join(print_columns)))
        
        for result in results:
            print_line = "{0}\t{1}\t{2}\t{3}".format(
                    result['ind_id'], result['fraction_of_common_variants'],
                    result['common_variants'], result['number_calls']
                )
            if outfile:
                outfile.write("{0}\n".format(print_line))
            else:
                print(print_line)
    
