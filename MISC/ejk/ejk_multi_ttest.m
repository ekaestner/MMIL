function ps = ejk_multi_ttest(struct1,struct2,cols)

ps = zeros(1,length(cols));

for i = 1:length(ps)
    [~,ps(i)] = ttest2(cell2mat(struct1(:,cols(i))),cell2mat(struct2(:,cols(i))));
end

end