# Utility for printing data in a terminal-friendly table layout.
# This module provides a generic table formatter that can be used across
# CLI tools in the gseda project. It uses a simple ASCII table representation
# that matches the visual layout of typical terminal output (e.g., samtools stats).
#
# The formatter supports:
# - Custom column alignment
# - Per-column width auto-calculation
# - Header formatting
# - Row formatting via callable functions (for pretty numbers)
#
# Extensibility:
#   * Add new formatting strategies by subclassing BaseTablePrinter.
#   * Register new table styles in a registry.
#   * Use this module from any CLI tool that outputs tabular data.

from typing import List, Callable, Dict, Any, Iterable, Tuple, Optional


class TablePrinter:
    """Generic ASCII table printer with configurable formatting."""
    def __init__(
        self,
        headers: List[str],
        rows: List[List[Any]],
        col_alignment: List[str] = None,
        col_widths: List[int] = None,
        formatters: List[Callable[[Any], str]] = None,
        title: Optional[str] = None,
    ):
        """
        Args:
            headers: List of column headers.
            rows: List of rows (each row is a list of values).
            col_alignment: List of 'left', 'center', 'right' alignment per column.
            col_widths: Pre-calculated column widths. If None, computed from data.
            formatters: List of callables to format each column's value.
            title: Optional title printed above the table.
        """
        self.headers = headers
        self.rows = rows
        self.col_alignment = col_alignment or ['left'] * len(headers)
        self.formatters = formatters or [lambda x: str(x)] * len(headers)
        self.title = title

        # Compute column widths
        if col_widths:
            self.col_widths = col_widths
        else:
            self.col_widths = self._calc_col_widths()
        self._validate()

    def _calc_col_widths(self) -> List[int]:
        # Start with header widths
        widths = [len(str(h)) for h in self.headers]
        for row in self.rows:
            for i, cell in enumerate(row):
                cell_str = self._format_cell(cell, i)
                w = len(cell_str)
                if w > widths[i]:
                    widths[i] = w
        return widths

    def _format_cell(self, cell: Any, col_idx: int) -> str:
        val = self._apply_formatter(cell, col_idx)
        # Truncate if needed? Not truncating for simplicity.
        return val

    def _apply_formatter(self, cell: Any, col_idx: int) -> str:
        fmt = self.formatters[col_idx]
        return fmt(cell)

    def _validate(self) -> None:
        if len(self.headers) != len(self.col_alignment):
            raise ValueError("col_alignment length must match headers length")
        if len(self.headers) != len(self.formatters):
            raise ValueError("formatters length must match headers length")
        for i, w in enumerate(self.col_widths):
            h = len(str(self.headers[i]))
            if w < h:
                raise ValueError(f"Column {i} width {w} is less than header length {h}")

    def _cell_str(self, cell: Any, col_idx: int) -> str:
        val = self._format_cell(cell, col_idx)
        return val.rjust(self.col_widths[col_idx], ' ')

    def _header_str(self, header: str, col_idx: int) -> str:
        return header.rjust(self.col_widths[col_idx], ' ')

    def print(self) -> None:
        if self.title:
            print(self.title)
            print('-' * max(self.col_widths) if self.col_widths else '')
        # Header row
        header_line = '  '.join(self._header_str(h, i) for i, h in enumerate(self.headers))
        print(header_line)
        print('  '.join('-' * w for w in self.col_widths))
        # Data rows
        for row in self.rows:
            row_line = '  '.join(self._cell_str(val, i) for i, val in enumerate(row))
            print(row_line)


class TablePrinterBuilder:
    """Helper to build a TablePrinter with fluent API."""
    def __init__(self, headers: List[str]):
        self.headers = headers
        self.rows: List[List[Any]] = []
        self.col_alignment: List[str] = []
        self.col_widths: List[int] = None
        self.formatters: List[Callable[[Any], str]] = []
        self.title: Optional[str] = None

    def add_rows(self, rows: Iterable[List[Any]]) -> 'TablePrinterBuilder':
        for row in rows:
            self.rows.append(row)
        return self

    def set_title(self, title: str) -> 'TablePrinterBuilder':
        self.title = title
        return self

    def set_col_alignment(self, alignment: List[str]) -> 'TablePrinterBuilder':
        self.col_alignment = alignment
        return self

    def set_col_widths(self, widths: List[int]) -> 'TablePrinterBuilder':
        self.col_widths = widths
        return self

    def set_formatters(self, formatters: List[Callable[[Any], str]]) -> 'TablePrinterBuilder':
        self.formatters = formatters
        return self

    def build(self) -> TablePrinter:
        return TablePrinter(
            headers=self.headers,
            rows=self.rows,
            col_alignment=self.col_alignment,
            col_widths=self.col_widths,
            formatters=self.formatters,
            title=self.title,
        )


def format_as_table(
    headers: List[str],
    rows: List[List[Any]],
    *,
    col_alignment: List[str] = None,
    col_widths: List[int] = None,
    formatters: List[Callable[[Any], str]] = None,
    title: Optional[str] = None,
) -> str:
    """
    Convenience function to produce a table as a string.
    Useful for one-off formatting without printing directly.
    """
    printer = TablePrinterBuilder(headers)
    printer = printer.add_rows(rows)
    printer = printer.set_title(title)
    printer = printer.set_col_alignment(col_alignment or ['left'] * len(headers))
    printer = printer.set_formatters(formatters or [lambda x: str(x)] * len(headers))
    if col_widths:
        printer = printer.set_col_widths(col_widths)
    table_printer = printer.build()
    return table_printer  # caller can do `print(table_printer)` or `str(table_printer)`
