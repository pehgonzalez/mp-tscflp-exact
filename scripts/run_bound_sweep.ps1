<#
.SYNOPSIS
    Measures the root bounds of the whole benchmark and survives being
    interrupted at any point, including a power loss.

.DESCRIPTION
    The study states the hierarchy v_BR = v_LP <= v_LD <= v* as a proposition
    and never measures the first three terms. This sweep measures them. It
    calls python\measure_bounds.py once per instance, and that script appends
    one row to logs\bounds.csv carrying

      v_lp      LP relaxation of the compact model
      v_br      LP value of the Benders reformulation, root cut loop run to
                convergence with the location variables continuous
      root_mip  dual bound at the end of the root node of the compact MIP
      root_bbc  dual bound at the end of the root node of the Benders master

    v_LD is not measured here. Every canonical L-BBC row of logs\results.csv
    already carries it in lag_lb, so recomputing it under a second
    implementation would put two numbers for one quantity in one paper.

    The sweep touches neither the solver binary nor logs\results.csv, so the
    campaign the paper reports stays exactly as it was produced.

    What the script gives beyond a plain loop over the instance folder.

      - it reports its own state on every launch, meaning how many rows are
        already recorded, how many are pending and how much machine time the
        remainder is expected to cost
      - it resumes from the point where it stopped, because the pending set is
        rebuilt from bounds.csv itself rather than from any private bookkeeping
      - measure_bounds.py builds its row in memory and appends it once, at the
        very end, so an interruption leaves bounds.csv holding whole rows only
      - it repairs a bounds.csv whose last line was truncated by a power cut
        before parsing it, moving the damaged text aside rather than deleting
      - it keeps a journal of every attempt, so an instance that starts and
        never finishes stays visible and is retried a bounded number of times
        instead of forever
      - it holds a lock carrying both the process id and that process's start
        time, so a lock left behind by a reboot is reclaimed automatically
        while a genuinely running sweep is never displaced
      - it enforces a wall-clock watchdog per instance, so a measurement that
        hangs past its own cap cannot stall the sweep
      - it stops cleanly between instances when a STOP file appears
      - it checks the interpreter and the licence before it starts, because a
        size-limited licence accepts the self test and then refuses every real
        instance, which would otherwise show up a hundred times in a row
      - it can register itself as a logon or startup task so an unattended
        machine resumes on its own after a reboot

    The estimate of the remaining time is calibrated on the sweep's own
    measured pace, per size class, and falls back to the worst case of one
    capped measurement per requested quantity while no pace has been observed.
    Running -Pilot first, which measures one instance of each size class,
    replaces that worst case with a real figure in about twenty minutes.

.PARAMETER Plan
    Print the pending list and the time estimate, run nothing.

.PARAMETER Status
    Print progress, the quality of what is already recorded, the lock state
    and any abandoned instances, run nothing.

.PARAMETER Pilot
    Measure one instance of each of the four size classes and stop. This is
    what turns the worst-case estimate into a calibrated one.

.PARAMETER Preflight
    Resolve the interpreter, check gurobipy and its licence, check that the
    two modules measure_bounds.py imports are importable, then stop.

.PARAMETER Selftest
    Run measure_bounds.py --selftest, which builds a tiny instance and checks
    that v_lp and v_br agree and that neither root bound falls below v_lp.
    Small enough to pass under a size-limited licence, so it verifies the
    measurement code and not the machine.

.PARAMETER Cap
    Wall-clock cap of each single measurement, in seconds. Four measurements
    per instance means the cap bounds one instance at four times this value
    plus model construction.

.PARAMETER Measures
    Comma-separated subset of v_lp, v_br, root_mip, root_bbc. Narrowing it
    both shortens the sweep and narrows what counts as a complete row.

.PARAMETER CorePoint
    Core-point cuts inside the measurement of root_bbc, off by default on
    every class. The table compares bound strength across classes, so a
    setting that varied with the class would confound the comparison the table
    exists to make. The value is recorded in the core_point column and forms
    part of the resume key, so a second sweep at the other setting adds rows
    rather than colliding with these.

.PARAMETER RerunIncomplete
    Also queue instances whose recorded row is missing a requested measure or
    carries a zero in its ok column. Use it after raising -Cap. bounds.csv is
    append-only, so the new row joins the old one and the analysis reads the
    later of the two. The per-instance attempt ceiling still applies, so this
    cannot loop.

.PARAMETER MaxRuns
    Stop after this many instances in the current session, leaving the rest
    for the next launch. Zero means no limit.

.PARAMETER Install
    Register a scheduled task that relaunches this script after a reboot.
    Uses ONLOGON, which needs no administrator rights. Add -AtStartup for a
    task that fires before logon, which does need them.

.PARAMETER Uninstall
    Remove that scheduled task.

.EXAMPLE
    .\run_bound_sweep.ps1 -Preflight
    .\run_bound_sweep.ps1 -Selftest
    .\run_bound_sweep.ps1 -Pilot
    .\run_bound_sweep.ps1 -Plan
    .\run_bound_sweep.ps1
    .\run_bound_sweep.ps1 -Status

.NOTES
    PowerShell 5.1 or later. Run it from anywhere. The script locates the tree
    that holds python\model_gurobipy.py and logs\, which in this project is
    the repository root, sets its working directory there and prints what it
    resolved before doing anything. The benchmark lives inside that tree, in
    data\instances, and a candidate instance folder is accepted only once it
    is seen to hold all hundred PSC files. Every choice can be overridden
    with -WorkDir, -DataDir, -Measure, -Out and -Python.

    To stop a sweep that is already running, create the file
    logs\boundsweep\STOP. The current instance finishes and the script exits
    between instances, leaving nothing half written.
#>

[CmdletBinding()]
param(
    [string]   $WorkDir,
    [string]   $DataDir,
    [string]   $Measure,
    [string]   $Out,
    [string]   $Python,
    [double]   $Cap            = 300.0,
    [int]      $Threads        = 16,
    [int]      $Seed           = 0,
    [int]      $CorePoint      = 0,
    [string]   $Measures       = 'v_lp,v_br,root_mip,root_bbc',
    [string]   $InstanceFilter = '*',
    [switch]   $RerunIncomplete,
    [int]      $MaxRuns        = 0,
    [int]      $MaxAttempts    = 2,
    [int]      $GraceSeconds   = 600,
    [switch]   $Plan,
    [switch]   $Status,
    [switch]   $Pilot,
    [switch]   $Preflight,
    [switch]   $Selftest,
    [switch]   $Install,
    [switch]   $Uninstall,
    [switch]   $AtStartup,
    [string]   $TaskName       = 'MP-TSCFLP-boundsweep'
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version 2.0

# ---------------------------------------------------------------- paths ----

$Here = $PSScriptRoot
if (-not $Here) { $Here = (Get-Location).Path }
$Parent = (Resolve-Path (Join-Path $Here '..')).Path

# The sweep runs in the tree that holds the python folder and the log folder,
# which in this project is the repository root. Rather than assume a layout,
# each candidate is tested for both.
$treeCands = @($Parent, $Here)
if (-not $WorkDir) {
    foreach ($c in $treeCands) {
        if ((Test-Path (Join-Path $c 'python\model_gurobipy.py')) -and
            (Test-Path (Join-Path $c 'logs'))) { $WorkDir = $c; break }
    }
}
if (-not $WorkDir) {
    throw ("No working tree found. None of the folders searched holds both " +
           "python\model_gurobipy.py and logs\, so the sweep would write its " +
           "bounds file somewhere the study does not read. Pass -WorkDir. " +
           "Searched: " + ($treeCands -join ' ; '))
}
$WorkDir = (Resolve-Path $WorkDir).Path

if (-not $Measure) { $Measure = Join-Path $WorkDir 'python\measure_bounds.py' }
if (-not (Test-Path $Measure)) {
    throw ("measure_bounds.py not found at $Measure. Copy it into " +
           (Join-Path $WorkDir 'python') + " or pass -Measure.")
}
$Measure = (Resolve-Path $Measure).Path
$PyDir   = Split-Path -Parent $Measure

# The hundred PSC instances live in data\instances since the 2026-08
# reorganisation, and a candidate is accepted only when it holds the whole
# benchmark.
$dataCands = @((Join-Path $WorkDir 'data\instances'),
               (Join-Path $WorkDir 'data'))
if (-not $DataDir) {
    foreach ($c in $dataCands) {
        if (-not (Test-Path $c)) { continue }
        $n = @(Get-ChildItem -Path $c -Filter 'PSC*.txt' -File -ErrorAction SilentlyContinue).Count
        if ($n -ge 100) { $DataDir = $c; break }
    }
}
if (-not $DataDir) {
    throw ("No instance folder found holding the hundred PSC*.txt files. " +
           "Pass -DataDir. Searched: " + ($dataCands -join ' ; '))
}
$DataDir = (Resolve-Path $DataDir).Path

$ScriptPath = $MyInvocation.MyCommand.Path
if (-not $ScriptPath) { $ScriptPath = Join-Path $Here 'run_bound_sweep.ps1' }

$LogDir   = Join-Path $WorkDir 'logs'
$StateDir = Join-Path $LogDir  'boundsweep'
if (-not $Out) { $Out = Join-Path $LogDir 'bounds.csv' }
$BoundsCsv   = $Out
$Journal     = Join-Path $StateDir 'journal.csv'
$LockFile    = Join-Path $StateDir 'sweep.lock'
$StopFile    = Join-Path $StateDir 'STOP'
$LauncherCmd = Join-Path $StateDir 'resume.cmd'
$PreflightPy = Join-Path $StateDir 'preflight.py'

# The header measure_bounds.py writes. Kept here so a bounds.csv from an older
# revision of that script is refused rather than silently half read.
$CsvHeader = 'datetime,instance,I,J,K,L,seed,threads,cap_s,core_point,' +
             'v_lp,v_lp_s,v_lp_ok,' +
             'v_br,v_br_s,v_br_ok,v_br_rounds,v_br_cuts,' +
             'root_mip,root_mip_s,root_mip_ok,' +
             'root_bbc,root_bbc_s,root_bbc_ok,root_bbc_cuts,' +
             'gurobi_version'
$CsvFieldCount = ($CsvHeader -split ',').Count

$Wanted = @($Measures -split ',' | ForEach-Object { $_.Trim() } | Where-Object { $_ -ne '' })
$Known  = @('v_lp', 'v_br', 'root_mip', 'root_bbc')
foreach ($w in $Wanted) {
    if ($Known -notcontains $w) {
        throw ("Unknown measure '$w'. Choose from " + ($Known -join ', ') + '.')
    }
}
if ($Wanted.Count -eq 0) { throw 'The -Measures list is empty, so there is nothing to measure.' }
$Measures = ($Wanted -join ',')

function Inv([object]$v) {
    # The cap leaves this script as text in three places, the command line of
    # the measurement, the journal and the launcher batch file, and all three
    # are read back by something that expects a dot. A bare interpolation of a
    # floating point number writes a comma wherever the machine culture uses
    # one, which is the same trap that made [double]::TryParse unusable here.
    if ($v -is [double] -or $v -is [single] -or $v -is [decimal]) {
        return ([double]$v).ToString('0.####', [System.Globalization.CultureInfo]::InvariantCulture)
    }
    return [string]$v
}
$CapText = Inv $Cap

New-Item -ItemType Directory -Force -Path $LogDir, $StateDir | Out-Null

function Say([string]$m)  { Write-Host ("[{0}] {1}" -f (Get-Date -Format 'HH:mm:ss'), $m) }
function Warn([string]$m) { Write-Warning $m }

function Note([string]$evt, [hashtable]$f) {
    # One journal line per state change. Append-only and closed after each
    # write, so a power cut loses at most the line being written.
    if (-not (Test-Path $Journal)) {
        Set-Content -Path $Journal -Encoding ASCII `
            -Value 'utc,event,instance,seed,threads,core_point,cap_s,attempt,detail'
    }
    $line = '{0},{1},{2},{3},{4},{5},{6},{7},{8}' -f `
        (Get-Date).ToUniversalTime().ToString('yyyy-MM-ddTHH:mm:ssZ'),
        $evt, $f['instance'], $Seed, $Threads, $CorePoint, $CapText, $f['attempt'],
        ([string]$f['detail'] -replace ',', ';')
    Add-Content -Path $Journal -Value $line -Encoding ASCII
}

# ------------------------------------------------------------ csv repair ----

function Assert-CsvHeader {
    # A header that does not match means the file was written by a different
    # revision of measure_bounds.py. Stopping is the only safe move, since the
    # repair below would otherwise treat every line of it as damaged.
    if (-not (Test-Path $BoundsCsv)) { return }
    $first = ''
    $reader = [System.IO.File]::OpenText($BoundsCsv)
    try { $first = [string]$reader.ReadLine() } finally { $reader.Close() }
    if ($null -eq $first) { return }
    if ($first.Trim() -ne $CsvHeader) {
        throw ("$BoundsCsv does not carry the expected $CsvFieldCount field header. " +
               "Move it aside before running the sweep, or point -Out elsewhere.")
    }
}

function Repair-BoundsCsv {
    # A power cut during the final append can leave a truncated row, which
    # would either break Import-Csv or parse into a record that makes the
    # sweep believe an instance had been measured. Drop any line that does not
    # carry the full field count, keeping a copy of everything.
    if (-not (Test-Path $BoundsCsv)) { return }
    $lines = [System.IO.File]::ReadAllLines($BoundsCsv)
    if ($lines.Count -eq 0) { return }

    $keep = New-Object System.Collections.Generic.List[string]
    $bad  = New-Object System.Collections.Generic.List[string]
    foreach ($l in $lines) {
        if ($l.Trim().Length -eq 0) { continue }
        if (($l -split ',').Count -eq $CsvFieldCount) { $keep.Add($l) } else { $bad.Add($l) }
    }
    if ($bad.Count -eq 0) { return }

    $stamp  = (Get-Date).ToString('yyyyMMdd-HHmmss')
    $rescue = Join-Path $StateDir ("bounds_damaged_{0}.txt" -f $stamp)
    $backup = Join-Path $StateDir ("bounds_before_repair_{0}.csv" -f $stamp)
    Set-Content -Path $rescue -Value $bad -Encoding ASCII
    Copy-Item -Path $BoundsCsv -Destination $backup -Force
    Set-Content -Path $BoundsCsv -Value $keep -Encoding ASCII

    Warn (("bounds.csv carried {0} malformed line(s), most likely an interrupted append. " +
           "They were moved to {1} and the original was copied to {2}.") -f $bad.Count, $rescue, $backup)
    Note 'repair' @{ instance = ''; attempt = 0; detail = "dropped $($bad.Count) malformed line(s)" }
}

function Backup-BoundsCsv {
    # One copy per session, before anything is appended.
    if (-not (Test-Path $BoundsCsv)) { return }
    $stamp = (Get-Date).ToString('yyyyMMdd-HHmmss')
    Copy-Item -Path $BoundsCsv -Destination (Join-Path $StateDir ("bounds_session_{0}.csv" -f $stamp)) -Force
}

# --------------------------------------------------------------- the plan ----

function ToNum([object]$v) {
    # A PowerShell cast parses with the invariant culture, unlike
    # [double]::TryParse, which follows the culture of the machine and would
    # reject 1801.23 wherever the comma is the decimal mark. The measurement
    # always writes numbers with a dot, so the invariant reading is correct.
    if ($null -eq $v) { return $null }
    $s = ([string]$v).Trim()
    if ($s -eq '') { return $null }
    try { return [double]$s } catch { return $null }
}

function Norm([object]$v) {
    $d = ToNum $v
    if ($null -eq $d) { return ([string]$v).Trim() }
    return ([long][Math]::Round($d)).ToString()
}

function Key([object]$inst, [object]$seed, [object]$threads, [object]$cp) {
    # The cap is deliberately absent from the key. A row measured under a
    # smaller cap is still a row about the same instance under the same
    # configuration, and raising the cap to sharpen it is what
    # -RerunIncomplete is for.
    '{0}|{1}|{2}|{3}' -f ([string]$inst).Trim(), (Norm $seed), (Norm $threads), (Norm $cp)
}

function Size-Class([string]$stem) {
    $p = $stem -split '-'
    if ($p.Count -ge 4) { return ($p[2] + '-' + $p[3]) }
    return 'other'
}

function Field([object]$row, [string]$name) {
    # Import-Csv gives one property per header column, but a file written by an
    # older revision could be short of one. Strict mode turns a missing
    # property into a terminating error, so it is read through the property bag.
    $p = $row.PSObject.Properties[$name]
    if ($null -eq $p) { return '' }
    if ($null -eq $p.Value) { return '' }
    return ([string]$p.Value).Trim()
}

function Row-Complete([object]$row) {
    # A row counts as complete when every requested quantity is present and
    # its own ok column says the measurement finished rather than hitting the
    # cap. This is what -RerunIncomplete tests.
    foreach ($w in $Wanted) {
        if ((Field $row $w) -eq '') { return $false }
        if ((Field $row ($w + '_ok')) -ne '1') { return $false }
    }
    return $true
}

function Get-Instances {
    $all = New-Object System.Collections.Generic.List[psobject]
    foreach ($p in 1..5) {
        foreach ($c in 1..5) {
            foreach ($sz in '50-5', '50-10', '100-5', '100-10') {
                $stem = "PSC$p-C$c-$sz"
                if ($stem -notlike $InstanceFilter) { continue }
                $all.Add([pscustomobject]@{ instance = $stem; cls = $sz })
            }
        }
    }
    return $all
}

function Get-Recorded {
    # Everything bounds.csv already holds, keyed the way the plan is keyed. The
    # later of two rows for the same key wins, because bounds.csv is
    # append-only and a rerun is meant to supersede what it repeats.
    $rows = @{}
    if (-not (Test-Path $BoundsCsv)) { return $rows }
    foreach ($r in (Import-Csv -Path $BoundsCsv)) {
        $k = Key (Field $r 'instance') (Field $r 'seed') (Field $r 'threads') (Field $r 'core_point')
        if ($rows.ContainsKey($k)) {
            if ((Field $r 'datetime') -lt (Field $rows[$k] 'datetime')) { continue }
        }
        $rows[$k] = $r
    }
    return $rows
}

function Get-Attempts {
    # How many times each instance has been started, from the journal. One that
    # keeps dying before it writes a row is abandoned rather than retried
    # forever, which is what keeps an unattended sweep moving.
    $a = @{}
    if (-not (Test-Path $Journal)) { return $a }
    foreach ($r in (Import-Csv -Path $Journal)) {
        if ((Field $r 'event') -ne 'start') { continue }
        $k = Key (Field $r 'instance') (Field $r 'seed') (Field $r 'threads') (Field $r 'core_point')
        if ($a.ContainsKey($k)) { $a[$k] = $a[$k] + 1 } else { $a[$k] = 1 }
    }
    return $a
}

function Get-Pace {
    # Observed wall time per size class, taken from the journal rather than
    # from the per-measurement seconds in bounds.csv, so that model building
    # and interpreter start-up are inside the figure the estimate uses.
    $sum = @{}
    $cnt = @{}
    $gotSum = 0.0
    $gotWorst = 0.0
    if (Test-Path $Journal) {
        foreach ($r in (Import-Csv -Path $Journal)) {
            if ((Field $r 'event') -ne 'done') { continue }
            $d = Field $r 'detail'
            # Written as a match rather than a negated one, so that the capture
            # group this reads is filled by the comparison that succeeded.
            if (-not ($d -match 'wall=([0-9.]+)')) { continue }
            $w = ToNum $Matches[1]
            if ($null -eq $w -or $w -le 0) { continue }
            $cl = Size-Class (Field $r 'instance')
            if ($sum.ContainsKey($cl)) { $sum[$cl] = $sum[$cl] + $w; $cnt[$cl] = $cnt[$cl] + 1 }
            else                       { $sum[$cl] = $w;             $cnt[$cl] = 1 }
            $c = ToNum (Field $r 'cap_s')
            if ($null -eq $c -or $c -le 0) { $c = $Cap }
            $gotSum   += $w
            $gotWorst += $c * $Wanted.Count
        }
    }
    $mean = @{}
    foreach ($cl in $sum.Keys) { $mean[$cl] = $sum[$cl] / $cnt[$cl] }
    $ratio = $null
    if ($gotWorst -gt 0) { $ratio = $gotSum / $gotWorst }
    return @{ mean = $mean; ratio = $ratio; n = $cnt }
}

function Estimate-Seconds($items, $pace) {
    # Cost of an instance is the mean already observed on its own size class.
    # Failing that it is the worst case, one full cap per requested quantity,
    # scaled by the ratio of observed to worst case over the classes that have
    # been seen. Failing that too it is the bare worst case, which is what the
    # first launch reports and what -Pilot exists to replace.
    $worst = $Cap * $Wanted.Count
    $s = 0.0
    foreach ($it in $items) {
        if ($null -eq $it) { continue }
        $cl = $it.cls
        if ($pace.mean.ContainsKey($cl))  { $s += $pace.mean[$cl] }
        elseif ($null -ne $pace.ratio)    { $s += $worst * $pace.ratio }
        else                              { $s += $worst }
    }
    return $s
}

function Show-Span([double]$sec) {
    $t = [TimeSpan]::FromSeconds([Math]::Round($sec))
    # Floor, not [int]: a PowerShell int cast rounds, so 2.85 hours would
    # print as 3h 51m, and every layer that reads the line back and prints
    # it again would inflate the estimate by another hour.
    return ('{0}h {1:00}m' -f [Math]::Floor($t.TotalHours), $t.Minutes)
}

# ------------------------------------------------------------------ lock ----

function Get-ProcStamp([int]$processId) {
    # Identity of a process is the pair of its id and its start time. After a
    # reboot the recorded id may well belong to something else, and the start
    # time is what tells the two apart.
    if ($processId -le 0) { return $null }
    try {
        $p = Get-Process -Id $processId -ErrorAction Stop
        return $p.StartTime.ToUniversalTime().ToString('yyyy-MM-ddTHH:mm:ss.fffZ')
    } catch { return $null }
}

function Acquire-Lock {
    if (Test-Path $LockFile) {
        $raw = ''
        try { $raw = [string](Get-Content -Path $LockFile -Raw -ErrorAction Stop) } catch { $raw = '' }
        $old = 0; $oldStamp = ''
        if ($raw -match 'pid=(\d+)')       { $old = [int]$Matches[1] }
        if ($raw -match 'procstart=(\S+)') { $oldStamp = $Matches[1] }
        $live = Get-ProcStamp $old
        if ($live -and ($oldStamp -eq '' -or $live -eq $oldStamp)) {
            Say ("A sweep is already running in process $old. Nothing to do here.")
            Say ("Create $StopFile to make it finish the current instance and exit.")
            return $false
        }
        Say "Stale lock from process $old, most likely a reboot or a kill. Reclaiming it."
        Note 'lock-reclaimed' @{ instance = ''; attempt = 0; detail = "stale pid $old" }
    }
    Set-Content -Path $LockFile -Encoding ASCII -Value @(
        "pid=$PID",
        "procstart=$(Get-ProcStamp $PID)",
        "host=$env:COMPUTERNAME",
        "started=$((Get-Date).ToUniversalTime().ToString('yyyy-MM-ddTHH:mm:ssZ'))"
    )
    return $true
}

function Release-Lock {
    if (-not (Test-Path $LockFile)) { return }
    $raw = ''
    try { $raw = [string](Get-Content -Path $LockFile -Raw -ErrorAction Stop) } catch { $raw = '' }
    if ($raw -match "pid=$PID\b") { Remove-Item -Path $LockFile -Force -ErrorAction SilentlyContinue }
}

# ------------------------------------------------------------ interpreter ----

$PreflightSource = @'
"""Preflight for the bound sweep, written by run_bound_sweep.ps1.

Three questions, answered in the order that makes a failure readable. Is there
an interpreter with gurobipy in it. Is its licence large enough for a real
instance, which the size-limited licence shipped with the pip wheel is not.
Are the two modules measure_bounds.py imports importable from the folder the
sweep will run it in.
"""
import sys

try:
    import gurobipy as gp
except Exception as exc:
    sys.stdout.write('NOGUROBI %s\n' % exc)
    raise SystemExit(2)

sys.stdout.write('PYOK %d.%d.%d\n' % gp.gurobi.version())
sys.stdout.write('PYEXE %s\n' % sys.executable)

# The smallest instance of the benchmark has fifty plants and five products and
# builds a model two orders of magnitude larger than this probe. A licence that
# refuses four thousand rows would therefore refuse every instance, and it
# would do so one instance at a time, a hundred times, after the sweep had
# already started.
N = 4000
try:
    env = gp.Env(empty=True)
    env.setParam('OutputFlag', 0)
    env.start()
    m = gp.Model('probe', env=env)
    x = m.addVars(N, lb=0.0, ub=1.0)
    for i in range(N):
        m.addConstr(x[i] <= 1.0)
    m.setObjective(gp.quicksum(x[i] for i in range(N)))
    m.optimize()
    sys.stdout.write('LICOK %d rows\n' % N)
except gp.GurobiError as exc:
    sys.stdout.write('LICSMALL %s\n' % exc)
    raise SystemExit(3)

if len(sys.argv) > 1:
    sys.path.insert(0, sys.argv[1])
try:
    import model_gurobipy
    import benders_gurobipy
except Exception as exc:
    sys.stdout.write('NOMODULES %s\n' % exc)
    raise SystemExit(4)
sys.stdout.write('MODOK\n')
'@

function Write-Preflight {
    Set-Content -Path $PreflightPy -Value $PreflightSource -Encoding ASCII
}

function Try-Python([string]$exe, [string[]]$pre) {
    # Returns the preflight output when the interpreter runs at all, and null
    # when the command does not exist. A command that exists and fails the
    # preflight returns its output, because the caller has to report why.
    if (-not (Get-Command $exe -ErrorAction SilentlyContinue)) { return $null }
    $argl = @()
    if ($pre) { $argl += $pre }
    $argl += @($PreflightPy, $PyDir)
    $text = ''
    $code = -1
    try {
        $text = (& $exe $argl 2>&1 | Out-String)
        $code = $LASTEXITCODE
    } catch {
        return $null
    }
    return @{ exe = $exe; pre = $pre; text = $text; code = $code }
}

function Resolve-Python {
    Write-Preflight
    $cands = New-Object System.Collections.Generic.List[psobject]
    if ($Python) { $cands.Add([pscustomobject]@{ exe = $Python; pre = @() }) }
    $cands.Add([pscustomobject]@{ exe = 'python';  pre = @() })
    $cands.Add([pscustomobject]@{ exe = 'python3'; pre = @() })
    $cands.Add([pscustomobject]@{ exe = 'py';      pre = @('-3') })

    $seen = @()
    foreach ($c in $cands) {
        $r = Try-Python $c.exe $c.pre
        if ($null -eq $r) { continue }
        $seen += $r
        if ($r.text -match 'MODOK') { return $r }
    }
    if ($seen.Count -eq 0) {
        throw ("No Python interpreter found. Tried python, python3 and py -3. " +
               "Pass -Python with the full path to the one that carries gurobipy.")
    }
    $first = $seen[0]
    $why = 'the preflight did not reach the end'
    if     ($first.text -match 'NOGUROBI')  { $why = 'gurobipy is not installed in it' }
    elseif ($first.text -match 'LICSMALL')  { $why = 'its Gurobi licence is size limited and refuses a four thousand row model, which every instance of the benchmark exceeds' }
    elseif ($first.text -match 'NOMODULES') { $why = "it cannot import model_gurobipy and benders_gurobipy from $PyDir" }
    throw ("No usable Python interpreter. The first one found was $($first.exe) and $why. " +
           "Full preflight output follows.`r`n" + $first.text)
}

# ------------------------------------------------------- scheduled task ----

function Write-Launcher {
    # schtasks caps the /TR string at 261 characters, and the full command line
    # with several quoted absolute paths is longer than that. A one line batch
    # file in the state folder keeps the registered command short. Every option
    # that shapes the plan is written out, not just the ones that differ from
    # the defaults, because a resume that silently reverted to a default would
    # sweep a different configuration from the one that was interrupted.
    $extra = ''
    if ($Python)          { $extra += (' -Python "{0}"' -f $Python) }
    if ($RerunIncomplete) { $extra += ' -RerunIncomplete' }
    $fmt = 'powershell.exe -NoProfile -ExecutionPolicy Bypass -WindowStyle Hidden -File "{0}" ' +
           '-WorkDir "{1}" -DataDir "{2}" -Measure "{3}" -Out "{4}" ' +
           '-Cap {5} -Threads {6} -Seed {7} -CorePoint {8} -Measures {9} ' +
           '-InstanceFilter "{10}" -MaxAttempts {11}{12}'
    $cmd = $fmt -f $ScriptPath, $WorkDir, $DataDir, $Measure, $BoundsCsv,
                   $CapText, $Threads, $Seed, $CorePoint, $Measures,
                   $InstanceFilter, $MaxAttempts, $extra
    $stamp = (Get-Date).ToString('yyyy-MM-dd HH:mm')
    $body = @(
        '@echo off',
        "rem Written by run_bound_sweep.ps1 -Install on $stamp",
        ('cd /d "{0}"' -f $WorkDir),
        $cmd
    )
    Set-Content -Path $LauncherCmd -Value $body -Encoding ASCII
    return $LauncherCmd
}

function Install-Task {
    $launcher = Write-Launcher
    $sched = 'ONLOGON'
    if ($AtStartup) { $sched = 'ONSTART' }
    Say "Registering scheduled task '$TaskName' with trigger $sched."
    Say "  launcher: $launcher"
    if ($AtStartup) {
        & schtasks.exe /Create /TN $TaskName /TR ('"' + $launcher + '"') /SC $sched /DELAY 0000:02 /RL HIGHEST /F
    } else {
        & schtasks.exe /Create /TN $TaskName /TR ('"' + $launcher + '"') /SC $sched /DELAY 0000:02 /F
    }
    if ($LASTEXITCODE -ne 0) {
        Warn ("schtasks returned $LASTEXITCODE. ONSTART needs an elevated prompt. " +
              "Re-run without -AtStartup for a logon trigger, which does not.")
    } else {
        Say "Done. The sweep resumes by itself after a reboot. Remove it with -Uninstall."
    }
}

function Uninstall-Task {
    & schtasks.exe /Delete /TN $TaskName /F
    if ($LASTEXITCODE -eq 0) { Say "Scheduled task '$TaskName' removed." }
}

# ----------------------------------------------------------- one instance ----

function Invoke-Measure($item, [int]$attempt, $py) {
    $inst = Join-Path $DataDir ($item.instance + '.txt')
    if (-not (Test-Path $inst)) { throw "Instance file not found: $inst" }

    $stamp = (Get-Date).ToString('yyyyMMdd-HHmmss')
    # These files stay in the state folder rather than in logs\, because the
    # provenance index of the campaign collects logs\console_*.txt and the
    # bound sweep is not a campaign run. Putting them anywhere under logs\ with
    # a name of that shape would add rows to a provenance file that describes
    # something else entirely.
    $base = Join-Path $StateDir ('bs_{0}_s{1}_cp{2}_{3}' -f $item.instance, $Seed, $CorePoint, $stamp)
    $out  = $base + '.out'
    $err  = $base + '.err'

    $argl = @()
    if ($py.pre) { $argl += $py.pre }
    # Only the three paths can carry a space, so only they are quoted. The rest
    # are numbers and short tags.
    $argl += @(('"{0}"' -f $Measure), ('"{0}"' -f $inst),
               '--cap', $CapText, '--threads', $Threads, '--seed', $Seed,
               '--papadakos', $CorePoint, '--measures', $Measures,
               '--out', ('"{0}"' -f $BoundsCsv))

    Note 'start' @{ instance = $item.instance; attempt = $attempt
                    detail = "measures=$Measures stamp=$stamp" }

    $t0 = Get-Date
    $p = Start-Process -FilePath $py.exe -ArgumentList $argl -WorkingDirectory $PyDir `
                       -NoNewWindow -PassThru `
                       -RedirectStandardOutput $out -RedirectStandardError $err
    # Touching the handle keeps ExitCode readable after the process is gone.
    $null = $p.Handle

    # Each measurement caps itself, so the sum of the caps bounds the solving.
    # The watchdog exists only for the case where one of them does not respect
    # its cap, and the grace window covers model construction, the subproblem
    # objects and interpreter start-up.
    $limitMs = [int]([Math]::Ceiling($Cap * $Wanted.Count) + $GraceSeconds) * 1000
    if (-not $p.WaitForExit($limitMs)) {
        Warn "Instance $($item.instance) exceeded its caps plus $GraceSeconds s. Terminating it."
        try { Stop-Process -Id $p.Id -Force -ErrorAction SilentlyContinue } catch { }
        Start-Sleep -Seconds 5
        Note 'killed' @{ instance = $item.instance; attempt = $attempt; detail = 'watchdog' }
        return $false
    }
    $wall = ((Get-Date) - $t0).TotalSeconds
    $code = -1
    try { $code = [int]$p.ExitCode } catch { $code = -1 }
    if ($code -ne 0) {
        Note 'exit-nonzero' @{ instance = $item.instance; attempt = $attempt
                               detail = "exit=$code wall=$([int]$wall)" }
        Warn "Instance $($item.instance) exited with code $code. See $err."
        return $false
    }

    # A zero exit with no row means the append itself failed. Treat it as a
    # failure so the instance is retried rather than silently lost.
    Repair-BoundsCsv
    $fresh = Get-Recorded
    $k = Key $item.instance $Seed $Threads $CorePoint
    if (-not $fresh.ContainsKey($k)) {
        Note 'no-row' @{ instance = $item.instance; attempt = $attempt
                         detail = "exit=0 but bounds.csv has no row wall=$([int]$wall)" }
        Warn "Instance $($item.instance) exited cleanly but wrote no row to bounds.csv."
        return $false
    }
    $full = Row-Complete $fresh[$k]
    $flag = 'complete'
    if (-not $full) { $flag = 'capped' }
    Note 'done' @{ instance = $item.instance; attempt = $attempt
                   detail = "wall=$([int]$wall) $flag" }
    if (-not $full) {
        Warn ("Instance $($item.instance) recorded a row in which at least one requested " +
              "measure hit the cap. Raise -Cap and relaunch with -RerunIncomplete to sharpen it.")
    }
    return $true
}

# ------------------------------------------------------------- entry ------

if ($Install)   { Install-Task;   exit 0 }
if ($Uninstall) { Uninstall-Task; exit 0 }

Set-Location $WorkDir
Assert-CsvHeader
Repair-BoundsCsv

if ($Preflight) {
    $py = Resolve-Python
    Say ("Interpreter  {0} {1}" -f $py.exe, ($py.pre -join ' '))
    foreach ($l in ($py.text -split "`r?`n")) { if ($l.Trim()) { Write-Host "    $l" } }
    Say 'Preflight passed. The interpreter, the licence and the two modules are all in place.'
    exit 0
}

if ($Selftest) {
    $py = Resolve-Python
    Say ("Interpreter  {0} {1}" -f $py.exe, ($py.pre -join ' '))
    $argl = @()
    if ($py.pre) { $argl += $py.pre }
    $argl += @(('"{0}"' -f $Measure), '--selftest')
    & $py.exe $argl
    $rc = $LASTEXITCODE
    if ($rc -eq 0) { Say 'Self test passed.' } else { Warn "Self test returned $rc." }
    exit $rc
}

$all      = @(Get-Instances)
$recorded = Get-Recorded
$attempts = Get-Attempts
$pace     = Get-Pace

$pending    = @()
$abandoned  = @()
$incomplete = @()
foreach ($it in $all) {
    $k = Key $it.instance $Seed $Threads $CorePoint
    if ($recorded.ContainsKey($k)) {
        if (Row-Complete $recorded[$k]) { continue }
        $incomplete += $it
        if (-not $RerunIncomplete) { continue }
    }
    $n = 0
    if ($attempts.ContainsKey($k)) { $n = $attempts[$k] }
    if ($n -ge $MaxAttempts) { $abandoned += $it } else { $pending += $it }
}
$pending    = @($pending)
$abandoned  = @($abandoned)
$incomplete = @($incomplete)

# What the sweep still owes, before -Pilot or -MaxRuns narrow the session. The
# progress line has to report this figure rather than the narrowed one, or a
# session capped at ten instances would claim ninety were already recorded.
$allPending = @($pending)

if ($Pilot) {
    # One instance of each size class, skipping any class whose pace is already
    # known, since a second measurement there would not sharpen the estimate.
    $pick = @()
    $seenCls = @()
    foreach ($it in $pending) {
        if ($seenCls -contains $it.cls) { continue }
        if ($pace.mean.ContainsKey($it.cls)) { continue }
        $seenCls += $it.cls
        $pick += $it
    }
    $pending = @($pick)
}
if ($MaxRuns -gt 0 -and $pending.Count -gt $MaxRuns) {
    $pending = @($pending[0..($MaxRuns - 1)])
}

$doneCount = $all.Count - $allPending.Count - $abandoned.Count
$eta = Estimate-Seconds $allPending $pace
$etaNow = Estimate-Seconds $pending $pace

Say ("Working tree {0}" -f $WorkDir)
Say ("Measurement  {0}" -f $Measure)
Say ("Instances    {0}" -f $DataDir)
Say ("Bounds file  {0}" -f $BoundsCsv)
Say ("Filter {0}   measures {1}   cap {2}s   threads {3}   seed {4}   core point {5}" -f `
     $InstanceFilter, $Measures, $CapText, $Threads, $Seed, $CorePoint)
Say ("Sweep: {0} of {1} instances recorded, {2} pending, {3} abandoned after {4} attempts." -f `
     $doneCount, $all.Count, $allPending.Count, $abandoned.Count, $MaxAttempts)
if ($incomplete.Count -gt 0) {
    $verb = 'queued for a rerun'
    if (-not $RerunIncomplete) { $verb = 'left alone, pass -RerunIncomplete to sharpen them' }
    Say ("{0} recorded row(s) hit the cap on at least one requested measure, {1}." -f `
         $incomplete.Count, $verb)
}
if ($pace.mean.Count -gt 0) {
    Say ("Estimated remaining machine time {0}, calibrated on {1} measured instance(s)." -f `
         (Show-Span $eta), (($pace.n.Values | Measure-Object -Sum).Sum))
} else {
    Say ("Estimated remaining machine time {0}. Nothing has been measured yet, so this is the " -f (Show-Span $eta))
    Say ("worst case of one full cap per quantity. Run -Pilot first to replace it with a real figure.")
}
if ($pending.Count -ne $allPending.Count) {
    Say ("This session is narrowed to {0} instance(s), about {1} of work." -f `
         $pending.Count, (Show-Span $etaNow))
}

if ($Status -or $Plan) {
    if ($Plan -and $pending.Count -gt 0) {
        Say 'Pending instances:'
        foreach ($q in $pending) { Write-Host ('    {0,-16} class {1}' -f $q.instance, $q.cls) }
    }
    if ($Status -and $recorded.Count -gt 0) {
        Say 'Quality of what is already recorded:'
        foreach ($w in $Wanted) {
            $okN = 0
            $noN = 0
            foreach ($r in $recorded.Values) {
                if ((Field $r ($w + '_ok')) -eq '1') { $okN++ } else { $noN++ }
            }
            Write-Host ('    {0,-9} finished {1,3}   hit the cap {2,3}' -f $w, $okN, $noN)
        }
    }
    if ($Status -and $pace.mean.Count -gt 0) {
        Say 'Measured pace per size class:'
        foreach ($cl in ($pace.mean.Keys | Sort-Object)) {
            Write-Host ('    {0,-8} {1,6:0} s per instance over {2} run(s)' -f `
                        $cl, $pace.mean[$cl], $pace.n[$cl])
        }
    }
    if ($abandoned.Count -gt 0) {
        Warn 'Abandoned. Inspect the .err files in logs\boundsweep before forcing a retry:'
        foreach ($q in $abandoned) { Write-Host ('    {0}' -f $q.instance) }
    }
    if ($Status -and (Test-Path $LockFile)) {
        Say 'Lock file present:'
        foreach ($l in (Get-Content $LockFile)) { Write-Host "    $l" }
    }
    exit 0
}

if ($pending.Count -eq 0) {
    if ($Pilot -and $allPending.Count -gt 0) {
        Say ('Every size class with pending work already has a measured pace, so the pilot ' +
             'has nothing to add. Relaunch without -Pilot to run the rest.')
    } else {
        Say 'Nothing left to measure.'
    }
    if ($abandoned.Count -gt 0) {
        Warn "$($abandoned.Count) instance(s) abandoned. Relaunch with a larger -MaxAttempts once the cause is fixed."
    }
    exit 0
}
if (Test-Path $StopFile) {
    Say "A STOP file is present at $StopFile. Delete it to let the sweep run."
    exit 0
}

$py = Resolve-Python
Say ("Interpreter  {0} {1}" -f $py.exe, ($py.pre -join ' '))

if (-not (Acquire-Lock)) { exit 0 }
Backup-BoundsCsv
Note 'sweep-start' @{ instance = ''; attempt = 0; detail = "pending=$($pending.Count)" }
$t0 = Get-Date
$ok = 0; $bad = 0; $i = 0

try {
    foreach ($it in $pending) {
        $i++
        if (Test-Path $StopFile) {
            Say 'STOP file found. Finishing here. Relaunch to continue from this point.'
            break
        }
        $rest = @()
        for ($j = $i - 1; $j -lt $pending.Count; $j++) { $rest += $pending[$j] }
        Say ('({0}/{1}) {2}   about {3} of work left' -f `
             $i, $pending.Count, $it.instance, (Show-Span (Estimate-Seconds $rest $pace)))

        $k = Key $it.instance $Seed $Threads $CorePoint
        $n = 0
        if ($attempts.ContainsKey($k)) { $n = $attempts[$k] }
        if (Invoke-Measure $it ($n + 1) $py) { $ok++ } else { $bad++ }
        # Refresh the pace as the sweep goes, so the countdown printed above
        # sharpens rather than repeating the figure of the first launch.
        $pace = Get-Pace
    }
} finally {
    Release-Lock
    $spent = ((Get-Date) - $t0).TotalSeconds
    Note 'sweep-stop' @{ instance = ''; attempt = 0
                         detail = "ok=$ok failed=$bad spent=$([int]$spent)s" }
    Say ('Session finished. {0} instance(s) measured, {1} failed, {2} of wall time spent.' -f `
         $ok, $bad, (Show-Span $spent))
    if ($bad -gt 0) { Say 'Relaunch to retry the failures, or use -Status to see what is left.' }
    if ($Pilot) { Say 'Pilot done. Relaunch with -Plan to see the calibrated estimate for the rest.' }
}
