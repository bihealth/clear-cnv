"""Tests for the misc module.

If external programs are called then there are two variants of the test.  One has the suffix
"_integration" and the actual external program is called (must have been installed, of course).
The other has the suffix "_unit" and calls to external programs are mocked out and only the
mocking results are tested.
"""

from unittest.mock import MagicMock
from types import SimpleNamespace

from clearCNV import misc


def test_which():
    assert misc.execs_available(["ls", "sort", "bedtools", "bedops"])
    assert not misc.execs_available(["invented-and-bogus"])


def test_merge_bedfile_integration(tmp_path):
    path_out = str(tmp_path / "out.bed")
    args = SimpleNamespace(infile="tests/testdata/test_merge_bed_input.bed", outfile=path_out)
    misc.merge_bedfile(args)
    with open("tests/testdata/test_merge_bed_output.bed", "rt") as f:
        expected = f.read()
    with open(path_out, "rt") as f:
        actual = f.read()
    assert expected == actual


def test_merge_bedfile_unit(mocker, tmp_path):
    mock_popen = MagicMock()
    mocker.patch("clearCNV.misc.subprocess.Popen", mock_popen)
    args = SimpleNamespace(
        infile="tests/testdata/test_merge_bed_input.bed", outfile=str(tmp_path / "out.bed")
    )
    misc.merge_bedfile(args)
    assert mock_popen.call_count == 2


def test_merge_bedfiles_integration(tmp_path):
    path_out = str(tmp_path / "out.bed")
    misc.merge_bedfiles(
        [
            "tests/testdata/test_merge_beds_input_0.bed",
            "tests/testdata/test_merge_beds_input_1.bed",
        ],
        path_out,
    )
    with open("tests/testdata/test_merge_beds_output.bed", "rt") as f:
        expected = f.read()
    with open(path_out, "rt") as f:
        actual = f.read()
    assert expected == actual


def test_merge_bedfiles_unit(mocker, tmp_path):
    mock_popen = MagicMock()
    mocker.patch("clearCNV.misc.subprocess.Popen", mock_popen)
    misc.merge_bedfiles(
        [
            "tests/testdata/test_merge_beds_input_0.bed",
            "tests/testdata/test_merge_beds_input_1.bed",
        ],
        str(tmp_path / "out.bed"),
    )
    assert mock_popen.call_count == 3
