A = randn(1000, 1000);  % Τυχαίος 1000x1000 πίνακας
x = randn(1000,1);       % Τυχαίο αρχικό διάνυσμα
v_init = randn(1000,1);  % Τυχαίο αρχικό ιδιοδιάνυσμα
N = 50;                  % Δοκίμασε με 50 ιδιοτιμές αρχικά

[lambdas, vectors] = myfun(A, N, x, v_init);
