% read reference file
refdata = readtable('automate_fba_data1.xlsx', 'Sheet', 'data', 'ReadRowNames', true);
% read model
model = xls2model('468rxn_model.xlsx');

% declare array for storing SGAM stoichiometry
sgamlist = [];

%iterate through each condition
for i = 1:size(refdata,1)

    %replace carb, protein, biomass biosynthesis reactions with references
    model = addReaction(model, 'Carbohydrate_syn', 'reactionFormula', refdata.Carbohydrate_syn{i});
    model = addReaction(model, 'Protein_syn', 'reactionFormula', refdata.Protein_syn{i});
    model = addReaction(model, 'BIO_R', 'reactionFormula', refdata.BIO_R{i});
    
    %change ub/lb
    model = changeRxnBounds(model, 'Ext1', refdata{i, 4}, 'b');

    %FBA
    tmp = optimizeCbModel(model, 'max');
    
    %parameter search using Bisection method
    %set initial value for x0 and x1
    x0 = 1;
    x1 = 50;
    
    %actual growth rate
    actual = refdata.BIO_drain_1_day_(i);

    % find until the difference between fitted GR and actual GR is less
    % than 0.0001
    
    while abs(tmp.f1 - actual) > 1e-4
        
        % change stoichiometric coef of SGAM to x0
        model0 = model;
        model0.S(94,351) = -x0;
        model0.S(88,351) = x0;
        model0.S(286,351) = x0;
        
        % change stoichiometric coef of SGAM to x1
        model1 = model;
        model1.S(94,351) = -x1;
        model1.S(88,351) = x1;
        model1.S(286,351) = x1;

        %GR of SGAM = x0
        b0 = optimizeCbModel(model0, 'max');
        %GR of SGAM = x1
        b1 = optimizeCbModel(model1, 'max');
        
        %from biomass ∝ 1/SGAM
        % if actual growth rate is located between b0 and b1 growth rate
        if (b0.f1 > actual) && (b1.f1 < actual)

            %find mean(x2) of x0 and x1 
            x2 = (x0 + x1)/2;
            
            % set SGAM to x2
            model.S(94,351) = -x2;
            model.S(88,351) = x2;
            model.S(286,351) = x2;

            %optimize x2
            tmp = optimizeCbModel(model, 'max');

            %if the growth rate from SGAM = x2 is less than actual
            if tmp.f1 < actual
                % set x1 = x2
                x1 = x2;
            else
                % if the growth rate from SGAM = x2 is greater than actual
                % set x0 = x2
                x0 = x2;
            end
        else
            % if b1 greater than actual, increase searching range by:
            if (b1.f1 > actual)
                % set x0 = x1
                x0 = x1;
                % multiply x1 by 2
                x1 = x1*2;
            else
                % if x0 is still too much
                % set x1 = x0
                x1 = x0;
                % set x0 = 0.01 (this value could be changed)
                x0 = 0.01;
            end
        end
    end

    % add SGAM value of condition "i" to sgamlist
    sgamlist(i,:) = model.S(88,351);
   
    % concatenate and save the results to "solutions.csv"
    if i == 1
        sol = tmp.x;
    else
        sol = cat(2, sol, tmp.x);
    end
    
end
writematrix(sol, 'solutions.csv')