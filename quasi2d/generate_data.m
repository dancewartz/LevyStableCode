clear
for saff_length = logspace(-3,5,9)
    simulate_2d
    save(['saff_length' num2str(saff_length) '.mat'],'u','params','ts');
end
