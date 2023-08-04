% adapted from: 
% https://se.mathworks.com/matlabcentral/answers/488735-verifyerror-in-a-script-based-unit-test
% [accessed: 2023-08-04]

function assertError(func, errID, varargin)

import matlab.unittest.diagnostics.Diagnostic;
import matlab.unittest.constraints.Throws;

% generate constraint met when given error ID is throwns
throws = Throws(errID);

% check that function throws the right error ID
passed = throws.satisfiedBy(func);
diagnostic_text = "";

% output own input as error message instead if wrong error is thrown
if ~passed
    diag = Diagnostic.join(varargin{:}, throws.getDiagnosticFor(func));
    arrayfun(@diagnose, diag);
    diagnostic_text = strjoin({diag.DiagnosticText},[newline newline]);
end

assert(passed, diagnostic_text); 