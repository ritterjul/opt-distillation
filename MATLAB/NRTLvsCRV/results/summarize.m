function [] = summarize(letter)
diff_all = zeros(10,1);
iter = 1;
for a = 1:3
    for b = a+1:4
        for c = b+1:5
            number = strcat(num2str(a),num2str(b),num2str(c));
            load(strcat(letter,'_',number),'diff')
            diff_all(iter) = diff;
            iter = iter + 1;
        end
    end
end
save(letter,'diff_all')
end