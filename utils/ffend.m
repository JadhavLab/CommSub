function sct = ffend(sct)

while ~isstruct(sct)
  sct=sct{end};
end