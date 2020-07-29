import csv
import itertools
import sys

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage

    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")

    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")

    #print("PEople : " , people)
    #print("probabilities" ,probabilities)
def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }

    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    parents={'father' :0, 'mother':0}

    prop,gene,trait,result=0,0,0,1
    for person in people:

        if people[person]['mother'] is None and people[person]['father'] is None:
            if person in two_genes:
                gene=2
                prop =PROBS['gene'][gene]
            elif person in one_gene:
                gene=1
                prop = PROBS['gene'][gene]
            elif person not in one_gene and person not in two_genes:
                gene=0
                prop = PROBS['gene'][gene]

        else:

            for p in parents:
                if people[person][p] in two_genes:
                    parents[p]=(1-PROBS['mutation'])
                elif people[person][p] in one_gene:
                    parents[p] = (1 - PROBS['mutation']) * 0.5
                else:
                    parents[p] = PROBS['mutation']
            if person in two_genes:
                gene=2
                prop = parents['father']*parents['mother']
            elif person in one_gene:
                gene=1
                prop = ((1-parents['father'])*parents['mother'])+((1-parents['mother'])*parents['father'])
            elif person not in one_gene and person not in two_genes:
                gene=0
                prop = (1-parents['father'])*(1-parents['mother'])
        if person in have_trait:
            trait =PROBS['trait'][gene][True]
        else:
            trait = PROBS['trait'][gene][False]
        result *=prop*trait





    return result



def update(probabilities, one_gene, two_genes, have_trait, p):
   for people in probabilities:
       if people in two_genes:
           probabilities[people]['gene'][2] +=p
       elif people in one_gene:
           probabilities[people]['gene'][1] +=p
       elif people not in one_gene and people not in two_genes:
           probabilities[people]['gene'][0] +=p
       if people in have_trait:
           probabilities[people]['trait'][True] +=p
       else:
           probabilities[people]['trait'][False] +=p



def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    for people in probabilities:

        Sum_g = sum(probabilities[people]['gene'].values())
        Sum_t = sum(probabilities[people]['trait'].values())
        for values in probabilities[people]['gene']:
            probabilities[people]['gene'][values] /=Sum_g
        for values in probabilities[people]['trait']:
            probabilities[people]['trait'][values] /=Sum_t



if __name__ == "__main__":
    main()
