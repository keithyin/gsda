import subprocess
import os
import pathlib
cur_filepath = os.path.abspath(__file__)  # noqa: E402
cur_filepath = pathlib.Path(cur_filepath)
project_path = cur_filepath.parent.parent.parent.parent

server_root = cur_filepath.parent.parent

project_path = str(project_path.absolute())
server_root_path = str(server_root.absolute())

print(project_path)


def run_claude_code(error_msg):
    prompt = f"""
Analyze this error and search the codebase:

{error_msg}

找到错误原因，如果是用户提交的表单有问题，请给出用户修改建议。如果是代码问题，请反馈给用户，让用户联系开发者。 (注意：提出的建议是给终端用户看的，语言要准确，指示要明确， 不要有废话。)

服务端的代码在 {server_root} .
前端的代码在 {server_root}/frontend/src .

注意：错误分析过程中，不要做任何代码修改操作。
"""
    print(prompt)

    env = os.environ.copy()

    # 👇 关键：在代码里写死
    env.update({
        "ANTHROPIC_AUTH_TOKEN" :"ollama",
        "ANTHROPIC_API_KEY": "",
        "ANTHROPIC_BASE_URL": "http://192.168.3.35:11434",
        
        "PATH": os.environ.get("PATH", "")
    })

    result = subprocess.run(
        [
            "claude",
            "--print",
            prompt
        ],
        capture_output=True,
        text=True,
        env=env,
        cwd=project_path
    )
    
    print("STDOUT:\n", result.stdout)
    print("STDERR:\n", result.stderr)
    print("RETURN CODE:", result.returncode)

    return result.stdout

if __name__ == "__main__":
    err_msg = """
                2026/04/23 02:22:59 - INFO - gsmm2-aligned-metric Version: 0.25.6
2026/04/23 02:22:59 - INFO - gsetl Version: 0.9.0
2026/04/23 02:22:59 - INFO - cmd: gsetl --outdir /path/to/file1-metric non-aligned-bam --bam /path/to/file1.bam -o /path/to/file1-metric/file1.non_aligned_fact.csv
exit: open bam error: /path/to/file1.bam
Traceback (most recent call last):
  File "<string>", line 1, in <module>
  File "/root/projects/gsda/gseda/src/gseda/ppl/reads_quality_stats_v3.py", line 773, in main_cli
    main(bam_file=bam, ref_fa=ref,
  File "/root/projects/gsda/gseda/src/gseda/ppl/reads_quality_stats_v3.py", line 707, in main
    non_aligned_metric_analysis(
  File "/root/projects/gsda/gseda/src/gseda/ppl/reads_quality_stats_v3.py", line 503, in non_aligned_metric_analysis
    df = pl.read_csv(
  File "/root/miniconda3/envs/server/lib/python3.8/site-packages/polars/_utils/deprecation.py", line 91, in wrapper
    return function(*args, **kwargs)
  File "/root/miniconda3/envs/server/lib/python3.8/site-packages/polars/_utils/deprecation.py", line 91, in wrapper
    return function(*args, **kwargs)
  File "/root/miniconda3/envs/server/lib/python3.8/site-packages/polars/_utils/deprecation.py", line 91, in wrapper
    return function(*args, **kwargs)
  File "/root/miniconda3/envs/server/lib/python3.8/site-packages/polars/io/csv/functions.py", line 500, in read_csv
    df = _read_csv_impl(
  File "/root/miniconda3/envs/server/lib/python3.8/site-packages/polars/io/csv/functions.py", line 646, in _read_csv_impl
    pydf = PyDataFrame.read_csv(
FileNotFoundError: No such file or directory (os error 2): /path/to/file1-metric/file1.non_aligned_fact.csv
    """
    
    res = run_claude_code(error_msg=err_msg)
    print(res)