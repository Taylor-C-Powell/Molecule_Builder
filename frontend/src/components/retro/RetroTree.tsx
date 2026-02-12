import { useEffect, useRef } from "react";
import * as d3 from "d3";
import type { RetroNodeResponse } from "@/api/types";
import { Card } from "@/components/ui/Card";

interface RetroTreeProps {
  tree: RetroNodeResponse;
}

interface TreeNode {
  smiles: string;
  is_purchasable: boolean;
  reaction?: string;
  children?: TreeNode[];
}

function toTreeData(node: RetroNodeResponse): TreeNode {
  return {
    smiles: node.smiles,
    is_purchasable: node.is_purchasable,
    reaction: node.best_disconnection?.reaction_name,
    children: node.children.length > 0 ? node.children.map(toTreeData) : undefined,
  };
}

function truncate(s: string, max: number): string {
  return s.length > max ? s.slice(0, max - 1) + "\u2026" : s;
}

export function RetroTree({ tree }: RetroTreeProps) {
  const svgRef = useRef<SVGSVGElement>(null);

  useEffect(() => {
    if (!svgRef.current) return;

    const data = toTreeData(tree);
    const root = d3.hierarchy(data);
    const width = 900;
    const nodeHeight = 60;
    const treeHeight = (root.height + 1) * nodeHeight * 2;

    const treeLayout = d3.tree<TreeNode>().size([width - 100, treeHeight - 60]);
    treeLayout(root);

    const svg = d3.select(svgRef.current);
    svg.selectAll("*").remove();
    svg.attr("viewBox", `0 0 ${width} ${treeHeight}`);

    const g = svg.append("g").attr("transform", "translate(50, 30)");

    // Links
    g.selectAll(".link")
      .data(root.links())
      .enter()
      .append("path")
      .attr("class", "link")
      .attr("fill", "none")
      .attr("stroke", "#404040")
      .attr("stroke-width", 1.5)
      .attr("d", d3.linkVertical<d3.HierarchyLink<TreeNode>, d3.HierarchyPointNode<TreeNode>>()
        .x((d) => d.x)
        .y((d) => d.y) as never);

    // Reaction labels on links
    g.selectAll(".reaction-label")
      .data(root.links().filter((l) => (l.source.data as TreeNode).reaction))
      .enter()
      .append("text")
      .attr("x", (d) => ((d.source as d3.HierarchyPointNode<TreeNode>).x + (d.target as d3.HierarchyPointNode<TreeNode>).x) / 2)
      .attr("y", (d) => ((d.source as d3.HierarchyPointNode<TreeNode>).y + (d.target as d3.HierarchyPointNode<TreeNode>).y) / 2 - 6)
      .attr("text-anchor", "middle")
      .attr("fill", "#737373")
      .attr("font-size", "9px")
      .text((d) => truncate((d.source.data as TreeNode).reaction ?? "", 20));

    // Nodes
    const nodes = g
      .selectAll(".node")
      .data(root.descendants())
      .enter()
      .append("g")
      .attr("transform", (d) => `translate(${d.x},${d.y})`);

    nodes
      .append("rect")
      .attr("x", -55)
      .attr("y", -14)
      .attr("width", 110)
      .attr("height", 28)
      .attr("rx", 6)
      .attr("fill", (d) => (d.data as TreeNode).is_purchasable ? "#052e16" : "#141414")
      .attr("stroke", (d) => (d.data as TreeNode).is_purchasable ? "#22c55e" : "#262626")
      .attr("stroke-width", 1);

    nodes
      .append("text")
      .attr("text-anchor", "middle")
      .attr("dy", "0.35em")
      .attr("fill", "#ededed")
      .attr("font-size", "10px")
      .attr("font-family", "monospace")
      .text((d) => truncate((d.data as TreeNode).smiles, 16));

  }, [tree]);

  return (
    <Card className="overflow-auto">
      <svg ref={svgRef} className="w-full min-h-[300px]" />
    </Card>
  );
}
