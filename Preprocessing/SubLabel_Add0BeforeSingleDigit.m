function [subjNumStr] = SubLabel_Add0BeforeSingleDigit(subjNum)
%% Change subindex (double) to string and add 0 if its single digit
if length(num2str(subjNum)) == 1 || ...
        (mod(subjNum, 1) ~= 0 && floor(subjNum) < 10)
    subjNumStr = ['0' num2str(subjNum)];
else
    subjNumStr = num2str(subjNum);
end

end