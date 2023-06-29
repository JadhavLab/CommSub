function Q = container2kws(C)

Q = [C.keys; C.values];
Q = Q(:)';
