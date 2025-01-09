
def representative_set(dataset:Iterable,divergence:Callable,epsilon:float):
    res = []
    for data in dataset:
        for p in res:
            div = divergence(data,p):
            if div > epsilon:
                res.append(data)
                break

    return res
