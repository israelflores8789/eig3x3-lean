module

import Tests.KnownAnswer
import Tests.Regression
import Tests.Certificates

/-!
# Tests.Main — test driver

Sequences the suites and reports an exit code: 0 if all suites pass, 1
otherwise. Wired up as the package's `testDriver`, so `lake test` builds
and runs this executable. A failing assertion aborts its suite (the
remaining suites still run) and fails the build.
-/

namespace Eig3x3.Tests

def runAll : IO UInt32 := do
  let suites : List (String × IO Unit) :=
    [("known-answer", runKnownAnswer),
     ("regression", runRegression),
     ("certificates", runCertificates)]
  let mut failed := 0
  for (name, suite) in suites do
    IO.println s!"== {name} =="
    try
      suite
    catch e =>
      failed := failed + 1
      IO.eprintln s!"FAIL {name}: {e}"
  if failed == 0 then
    IO.println "all test suites passed"
    return 0
  else
    IO.eprintln s!"{failed} of {suites.length} test suites failed"
    return 1

end Eig3x3.Tests

def main : IO UInt32 := Eig3x3.Tests.runAll
