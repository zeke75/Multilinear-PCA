function d=distance(a,b,g,k)

% Distance function (Table II in the paper), used in Matching function

switch k
    case 1
        d=norm(a-b,1);
    case 2
        d=norm(a-b,2);
    case 3
        d=-dot(a,b)/(sqrt(sum(a.^2))*sqrt(sum(b.^2)));
    case 4
        d=-dot(a./g,b);
    case 5
        d=norm((a-b)./g,1);
    case 6
        d=norm((a-b)./sqrt(g),2);
end
        