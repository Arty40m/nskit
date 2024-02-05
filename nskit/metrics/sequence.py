from ..algo import levdist



def levsim(x, y):
    return 1 - levdist(x, y)/max(len(x), len(y))
    

def sublevsim(x, y):
    return 1 - (levdist(x, y) - abs(len(x)-len(y)))/min(len(x), len(y))