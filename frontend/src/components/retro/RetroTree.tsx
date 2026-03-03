import { useCallback, useEffect, useRef, useState } from "react";
import * as d3 from "d3";
import type { RetroNodeResponse, DisconnectionResponse } from "@/api/types";
import { Card } from "@/components/ui/Card";
import { NodeDetail } from "@/components/retro/NodeDetail";

interface RetroTreeProps {
  tree: RetroNodeResponse;
}

interface TreeNode {
  smiles: string;
  is_purchasable: boolean;
  functional_groups: string[];
  reaction?: string;
  score?: number;
  disconnections?: DisconnectionResponse[];
  children?: TreeNode[];
  _children?: TreeNode[];
}

function toTreeData(node: RetroNodeResponse): TreeNode {
  return {
    smiles: node.smiles,
    is_purchasable: node.is_purchasable,
    functional_groups: node.functional_groups,
    reaction: node.best_disconnection?.reaction_name,
    score: node.best_disconnection?.score,
    disconnections: node.disconnections,
    children: node.children.length > 0 ? node.children.map(toTreeData) : undefined,
  };
}

function truncate(s: string, max: number): string {
  return s.length > max ? s.slice(0, max - 1) + "\u2026" : s;
}

function scoreColor(score: number): string {
  if (score >= 70) return "#22c55e";
  if (score >= 40) return "#eab308";
  return "#ef4444";
}

type D3HNode = d3.HierarchyPointNode<TreeNode>;

export function RetroTree({ tree }: RetroTreeProps) {
  const svgRef = useRef<SVGSVGElement>(null);
  const [selectedNode, setSelectedNode] = useState<TreeNode | null>(null);

  const render = useCallback(() => {
    if (!svgRef.current) return;

    const data = toTreeData(tree);
    const width = 960;
    const nodeW = 130;
    const nodeH = 32;
    const root = d3.hierarchy(data);

    // Count visible leaves for height calc
    const leafCount = root.leaves().length;
    const treeHeight = Math.max(400, leafCount * 56);

    const treeLayout = d3.tree<TreeNode>().size([width - 160, treeHeight - 80]);
    treeLayout(root);

    const svg = d3.select(svgRef.current);
    svg.selectAll("*").remove();
    svg
      .attr("width", width)
      .attr("height", treeHeight + 60)
      .attr("viewBox", `0 0 ${width} ${treeHeight + 60}`);

    // Zoom container
    const container = svg.append("g");
    const zoomBehavior = d3.zoom<SVGSVGElement, unknown>()
      .scaleExtent([0.3, 3])
      .on("zoom", (event: d3.D3ZoomEvent<SVGSVGElement, unknown>) => {
        container.attr("transform", event.transform.toString());
      });
    svg.call(zoomBehavior);
    // Initial offset
    container.attr("transform", "translate(80, 40)");
    svg.call(zoomBehavior.transform, d3.zoomIdentity.translate(80, 40));

    const duration = 300;

    // Links with animation
    const links = container
      .selectAll<SVGPathElement, d3.HierarchyLink<TreeNode>>(".link")
      .data(root.links())
      .enter()
      .append("path")
      .attr("class", "link")
      .attr("fill", "none")
      .attr("stroke", "#404040")
      .attr("stroke-width", 1.5)
      .attr("d", d3.linkVertical<d3.HierarchyLink<TreeNode>, D3HNode>()
        .x((d) => d.x)
        .y((d) => d.y) as never)
      .attr("opacity", 0);
    links.transition().duration(duration).attr("opacity", 1);

    // Score badges on links
    container
      .selectAll<SVGTextElement, d3.HierarchyLink<TreeNode>>(".score-label")
      .data(root.links().filter((l) => (l.source.data as TreeNode).score != null))
      .enter()
      .append("text")
      .attr("x", (d) => ((d.source as D3HNode).x + (d.target as D3HNode).x) / 2)
      .attr("y", (d) => ((d.source as D3HNode).y + (d.target as D3HNode).y) / 2 - 8)
      .attr("text-anchor", "middle")
      .attr("fill", (d) => scoreColor((d.source.data as TreeNode).score ?? 0))
      .attr("font-size", "9px")
      .attr("font-weight", "600")
      .text((d) => {
        const s = (d.source.data as TreeNode).score;
        return s != null ? s.toFixed(0) : "";
      });

    // Reaction labels on links
    container
      .selectAll<SVGTextElement, d3.HierarchyLink<TreeNode>>(".reaction-label")
      .data(root.links().filter((l) => (l.source.data as TreeNode).reaction))
      .enter()
      .append("text")
      .attr("x", (d) => ((d.source as D3HNode).x + (d.target as D3HNode).x) / 2)
      .attr("y", (d) => ((d.source as D3HNode).y + (d.target as D3HNode).y) / 2 + 4)
      .attr("text-anchor", "middle")
      .attr("fill", "#737373")
      .attr("font-size", "8px")
      .text((d) => truncate((d.source.data as TreeNode).reaction ?? "", 24));

    // Nodes
    const nodes = container
      .selectAll<SVGGElement, D3HNode>(".node")
      .data(root.descendants())
      .enter()
      .append("g")
      .attr("class", "node")
      .attr("transform", (d) => `translate(${d.x},${d.y})`)
      .style("cursor", "pointer")
      .on("click", (_event, d) => {
        setSelectedNode(d.data);
      });

    // Node rectangles -- wider, with color coding
    const halfW = nodeW / 2;
    const halfH = nodeH / 2;

    nodes
      .append("rect")
      .attr("x", -halfW)
      .attr("y", -halfH)
      .attr("width", nodeW)
      .attr("height", nodeH)
      .attr("rx", 6)
      .attr("fill", (d) => {
        const n = d.data as TreeNode;
        if (n.is_purchasable) return "#052e16";
        return "#141414";
      })
      .attr("stroke", (d) => {
        const n = d.data as TreeNode;
        if (n.is_purchasable) return "#22c55e";
        // Dead end: not purchasable and no children
        if (!n.children && !n._children) return "#6b7280";
        return "#262626";
      })
      .attr("stroke-width", (d) => {
        const n = d.data as TreeNode;
        if (!n.is_purchasable && !n.children && !n._children) return 1;
        return 1.5;
      })
      .attr("stroke-dasharray", (d) => {
        const n = d.data as TreeNode;
        if (!n.is_purchasable && !n.children && !n._children) return "4,2";
        return "none";
      });

    // Purchasable icon (small circle)
    nodes
      .filter((d) => (d.data as TreeNode).is_purchasable)
      .append("circle")
      .attr("cx", -halfW + 10)
      .attr("cy", 0)
      .attr("r", 3)
      .attr("fill", "#22c55e");

    // SMILES text
    nodes
      .append("text")
      .attr("text-anchor", "middle")
      .attr("x", (d) => (d.data as TreeNode).is_purchasable ? 5 : 0)
      .attr("dy", "0.35em")
      .attr("fill", "#ededed")
      .attr("font-size", "9px")
      .attr("font-family", "monospace")
      .text((d) => truncate((d.data as TreeNode).smiles, 18));

    // Collapse/expand indicator for nodes with children
    nodes
      .filter((d) => {
        const n = d.data as TreeNode;
        return !!(n.children || n._children);
      })
      .append("text")
      .attr("x", halfW - 10)
      .attr("y", 0)
      .attr("dy", "0.35em")
      .attr("text-anchor", "middle")
      .attr("fill", "#737373")
      .attr("font-size", "10px")
      .text((d) => (d.data as TreeNode).children ? "-" : "+");
  }, [tree]);

  useEffect(() => {
    render();
  }, [render]);

  return (
    <Card className="relative overflow-hidden">
      <svg ref={svgRef} className="w-full min-h-[400px]" />
      {selectedNode && (
        <NodeDetail node={selectedNode} onClose={() => setSelectedNode(null)} />
      )}
    </Card>
  );
}
