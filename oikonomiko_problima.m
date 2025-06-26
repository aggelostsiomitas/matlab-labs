clear;
clc;
A=[0.1,0.2,0.2,0.6,0.2;0.1,0.1,0.1,0.1,0.2;0.1,0.3,0.4,0.1,0.2;0.3,0.3,0.1,0.1,0.2;0.4,0.1,0.2,0.1,0.2];
% Κόστος της διαφήμισης (ετησίως)
advertising_cost = 40e6;

% Συνολικές πωλήσεις καφέ σε κιλά
total_sales_kg = 60e6;

% Κέρδος ανά κιλό καφέ
profit_per_kg = 1;

%αρχικη πρωτη στηλη Α
initial_distribution = [0.1; 0.1; 0.1; 0.3; 0.4];


%υπολογισμος αρχικων πωλησεων
initial_sales_distribution = A * initial_distribution;

% Υπολογισμός των πωλήσεων σε κιλά για κάθε μάρκα
initial_sales_kg = initial_sales_distribution * total_sales_kg;

% Υπολογισμός του κέρδους πριν την αλλαγή
initial_profit = sum(initial_sales_kg) * profit_per_kg; 

%υπολογισμος νεου κερδους με χρηση της νεας αλλαγης του πρακτορειου
%διαφημισεων
A_new = A;
A_new(:, 1) = [0.3; 0.1; 0.1; 0.2; 0.3]; 

new_sales_distribution = A_new * initial_distribution;
new_sales_kg = new_sales_distribution * total_sales_kg;
new_profit = sum(new_sales_kg) * profit_per_kg; 
profit_increase = new_profit - initial_profit;

disp('Αρχικά Κέρδη:');
disp(initial_profit);
disp('Νέα Κέρδη μετά την παρέμβαση διαφήμισης:');
disp(new_profit);
disp('Αύξηση Κερδών λόγω διαφήμισης:');
disp(profit_increase);

if profit_increase > advertising_cost
    disp('Ο κατασκευαστής της μάρκας B1 πρέπει να νοικιάσει τις υπηρεσίες του πρακτορείου διαφημίσεων.');
else
    disp('Ο κατασκευαστής της μάρκας B1 δεν πρέπει να νοικιάσει τις υπηρεσίες του πρακτορείου διαφημίσεων.');
end