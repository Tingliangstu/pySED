clear; close all;
tic;

%% set up some parameters related to the dump file

Total_atoms = 42240;        % number of atoms (have to  change the number of atoms in your system )

len_header=9;
N=Total_atoms+len_header;              % length of block

total_time_step = 1;                   % total number of dump point
in_file = 'dump.optimize';             %  dump file name
save_file_name = 'basis.in';           %  output  file
dump_columns = 6;

%% deal with these velocity data

atom_all = zeros(Total_atoms, dump_columns, total_time_step);

fileID = fopen(in_file);
formatSpec = '%n %n %n %n %n %n';

for i = 1:total_time_step

    C = textscan(fileID, formatSpec, N, 'HeaderLines', len_header, 'Delimiter', '\t');

    for  j = 1 : dump_columns

        atom_all(:,j,i)=C{1, j};                                                          % Just extract the velocity information for x, y, and z

    end
end

fclose(fileID);

%% output basis.in
unitcell_atoms = 88;
unitcell_number = Total_atoms/unitcell_atoms;  %% must be integer
 
fid = fopen(save_file_name, 'w');
fprintf(fid, 'Create basis.in for SED method by using LT_Codes\n');
fprintf(fid, 'atoms_ids unitcell_index basis_index mass_types\n');
atom_index = 1;

for n = 1 : unitcell_number

    for basis_index = 1 : unitcell_atoms

        fprintf(fid, '%d %d %d %g \n', atom_all(atom_index, 1), n, ...
                                                basis_index, atom_all(atom_index, 6));

        atom_index = atom_index + 1;
        
    end
end

fclose(fid);

toc;



